from datetime import datetime, timedelta
import multiprocessing as mp
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import (
    EarthLocation,
    get_body,
    solar_system_ephemeris,
)
from astropy.time import Time, TimeDelta
from astropy.utils import iers

from constants import (
    CELESTIAL_BODIES,
    START_TIME_TT,
    TRAJECTORY_TABLES_DIRNAME,
    WORKER_COUNT,
)

input_dir = Path(TRAJECTORY_TABLES_DIRNAME)
output_tsv_path = Path("moon_occultation_events.tsv")
SIMULATION_START = datetime.strptime(START_TIME_TT, "%Y-%m-%d %H:%M:%S")

# Adaptive refinement settings.
COARSE_TRIGGER_MOON_RADII = 2.0
COARSE_GAP_FACTOR = 1.5
REFINEMENT_PADDING_MINUTES = 60.0
REFINEMENT_STEP_MINUTES = 1.0
REFINEMENT_NOISE_STD_DEG = 2.0 / 60.0

iers.conf.iers_degraded_accuracy = "ignore"
solar_system_ephemeris.set("de440")

OBSERVATORY = EarthLocation(lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m)
START_TIME_ASTROPY = Time(START_TIME_TT, scale="tt")

BODY_METADATA = {
    body: {
        "x_limit": x_limit,
        "orbital_period": orbital_period,
        "apparent_radius_deg": apparent_radius_deg,
    }
    for body, x_limit, orbital_period, apparent_radius_deg in CELESTIAL_BODIES
}

HOUR_IN_DAYS = 1.0 / 24.0
MINUTE_IN_DAYS = 1.0 / 1440.0
MERGE_EPSILON_DAYS = 1.0 / 86400.0


def merge_time_windows(windows, gap_days):
    if not windows:
        return []

    sorted_windows = sorted(windows, key=lambda window: window[0])
    merged = [sorted_windows[0]]

    for start, end in sorted_windows[1:]:
        last_start, last_end = merged[-1]
        if start - last_end <= gap_days + MERGE_EPSILON_DAYS:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    return merged


def extend_window(window, padding_minutes):
    padding_days = padding_minutes * MINUTE_IN_DAYS
    return window[0] - padding_days, window[1] + padding_days


def interpolate_crossing_time(day_a, day_b, value_a, value_b, threshold):
    delta = value_b - value_a
    if delta == 0:
        return day_b
    ratio = (threshold - value_a) / delta
    ratio = np.clip(ratio, 0.0, 1.0)
    return day_a + (day_b - day_a) * ratio


def extract_contact_times(
    day_offsets, separations_deg, outer_threshold_deg, inner_threshold_deg
):
    external_ingress = None
    internal_ingress = None
    internal_egress = None
    external_egress = None

    for index in range(len(separations_deg) - 1):
        sep_a = separations_deg[index]
        sep_b = separations_deg[index + 1]
        day_a = day_offsets[index]
        day_b = day_offsets[index + 1]

        if external_ingress is None and sep_a > outer_threshold_deg >= sep_b:
            external_ingress = interpolate_crossing_time(
                day_a,
                day_b,
                sep_a,
                sep_b,
                outer_threshold_deg,
            )

        if internal_ingress is None and sep_a > inner_threshold_deg >= sep_b:
            internal_ingress = interpolate_crossing_time(
                day_a,
                day_b,
                sep_a,
                sep_b,
                inner_threshold_deg,
            )

        if internal_ingress is not None and internal_egress is None:
            if sep_a <= inner_threshold_deg < sep_b:
                internal_egress = interpolate_crossing_time(
                    day_a,
                    day_b,
                    sep_a,
                    sep_b,
                    inner_threshold_deg,
                )

        if external_ingress is not None and external_egress is None:
            if sep_a <= outer_threshold_deg < sep_b:
                external_egress = interpolate_crossing_time(
                    day_a,
                    day_b,
                    sep_a,
                    sep_b,
                    outer_threshold_deg,
                )

    return (
        external_ingress,
        internal_ingress,
        internal_egress,
        external_egress,
    )


def load_trajectory(body):
    trajectory_hdf_path = input_dir / f"{body}_trajectory.h5"
    if not trajectory_hdf_path.exists():
        raise FileNotFoundError(f"Missing input HDF5 for {body}: {trajectory_hdf_path}")

    trajectory_data = pd.read_hdf(trajectory_hdf_path, key="trajectory")
    if trajectory_data.empty:
        raise ValueError(f"Input HDF5 for {body} is empty: {trajectory_hdf_path}")

    required_columns = {"time_point", "ra_filtered", "dec_filtered"}
    missing = required_columns.difference(trajectory_data.columns)
    if missing:
        raise ValueError(f"Input HDF5 for {body} is missing columns: {sorted(missing)}")

    return trajectory_data[["time_point", "ra_filtered", "dec_filtered"]].copy()


def angular_separation_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    ra1_rad = np.radians(ra1_deg)
    dec1_rad = np.radians(dec1_deg)
    ra2_rad = np.radians(ra2_deg)
    dec2_rad = np.radians(dec2_deg)

    cos_sep = np.sin(dec1_rad) * np.sin(dec2_rad) + np.cos(dec1_rad) * np.cos(
        dec2_rad
    ) * np.cos(ra1_rad - ra2_rad)
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.degrees(np.arccos(cos_sep))


def wrapped_delta_ra_deg(ra_a_deg, ra_b_deg):
    return ((ra_a_deg - ra_b_deg + 180.0) % 360.0) - 180.0


def vector_angle_deg(vector_a, vector_b):
    norm_a = np.linalg.norm(vector_a)
    norm_b = np.linalg.norm(vector_b)
    if norm_a == 0.0 or norm_b == 0.0:
        return None

    cosine = np.dot(vector_a, vector_b) / (norm_a * norm_b)
    cosine = np.clip(cosine, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def contact_snapshot(
    contact_points,
    fine_days,
    contact_time,
):
    if contact_time is None:
        return None

    contact_index = int(np.argmin(np.abs(fine_days - contact_time)))
    return contact_points.iloc[contact_index]


def center_to_body_vector(snapshot):
    if snapshot is None:
        return None

    return np.array(
        [
            wrapped_delta_ra_deg(snapshot["body_ra_noisy"], snapshot["moon_ra_noisy"]),
            snapshot["body_dec_noisy"] - snapshot["moon_dec_noisy"],
        ],
        dtype=float,
    )


def day_offset_to_date_string(day_offset):
    total_seconds = int(np.round(float(day_offset) * 86400.0))
    timestamp = SIMULATION_START + timedelta(seconds=total_seconds)
    return timestamp.strftime("%Y-%m-%d %H:%M:%S")


def format_occultation_duration(start_day_offset, end_day_offset):
    if start_day_offset is None or end_day_offset is None:
        return ""

    duration_seconds = int(
        np.round((float(end_day_offset) - float(start_day_offset)) * 86400.0)
    )
    if duration_seconds < 0:
        return ""

    return str(duration_seconds)


def process_body(
    body,
    body_apparent_radius_deg,
    moon_data,
    moon_apparent_radius_deg,
):
    rng = np.random.default_rng(42 + sum(ord(char) for char in body))

    body_data = load_trajectory(body)
    coarse_data = body_data.merge(
        moon_data,
        on="time_point",
        how="inner",
        suffixes=("_body", "_moon"),
    )
    if coarse_data.empty:
        return []

    coarse_separation_deg = angular_separation_deg(
        coarse_data["ra_filtered_body"].to_numpy(),
        coarse_data["dec_filtered_body"].to_numpy(),
        coarse_data["ra_filtered_moon"].to_numpy(),
        coarse_data["dec_filtered_moon"].to_numpy(),
    )

    coarse_threshold_deg = COARSE_TRIGGER_MOON_RADII * moon_apparent_radius_deg
    coarse_hits = np.asarray(
        coarse_data.loc[
            coarse_separation_deg <= coarse_threshold_deg,
            "time_point",
        ],
        dtype=float,
    )
    if coarse_hits.size == 0:
        return []

    point_windows = [(float(day), float(day)) for day in np.sort(coarse_hits)]
    merged_windows = merge_time_windows(point_windows, gap_days=HOUR_IN_DAYS)

    body_radius_deg = body_apparent_radius_deg
    outer_threshold_deg = moon_apparent_radius_deg + body_radius_deg
    inner_threshold_deg = abs(moon_apparent_radius_deg - body_radius_deg)
    complete_occultation_possible = moon_apparent_radius_deg > body_radius_deg

    event_rows = []
    event_index = 0
    fine_step_days = REFINEMENT_STEP_MINUTES * MINUTE_IN_DAYS

    for coarse_start, coarse_end in merged_windows:
        fine_start, fine_end = extend_window(
            (coarse_start, coarse_end),
            REFINEMENT_PADDING_MINUTES,
        )

        fine_days = np.arange(
            fine_start,
            fine_end + fine_step_days / 2.0,
            fine_step_days,
        )
        fine_times = START_TIME_ASTROPY + TimeDelta(fine_days * u.day)

        moon_coords = get_body("moon", fine_times, location=OBSERVATORY)
        body_coords = get_body(body, fine_times, location=OBSERVATORY)

        moon_ra = np.asarray(moon_coords.ra.deg, dtype=float)
        moon_dec = np.asarray(moon_coords.dec.deg, dtype=float)
        body_ra = np.asarray(body_coords.ra.deg, dtype=float)
        body_dec = np.asarray(body_coords.dec.deg, dtype=float)

        moon_ra_noisy = moon_ra + rng.normal(
            0.0,
            REFINEMENT_NOISE_STD_DEG,
            size=len(fine_days),
        )
        moon_dec_noisy = moon_dec + rng.normal(
            0.0,
            REFINEMENT_NOISE_STD_DEG,
            size=len(fine_days),
        )
        body_ra_noisy = body_ra + rng.normal(
            0.0,
            REFINEMENT_NOISE_STD_DEG,
            size=len(fine_days),
        )
        body_dec_noisy = body_dec + rng.normal(
            0.0,
            REFINEMENT_NOISE_STD_DEG,
            size=len(fine_days),
        )

        contact_points = pd.DataFrame(
            {
                "body_ra_noisy": body_ra_noisy,
                "body_dec_noisy": body_dec_noisy,
                "moon_ra_noisy": moon_ra_noisy,
                "moon_dec_noisy": moon_dec_noisy,
            }
        )

        fine_separation_deg = angular_separation_deg(
            body_ra_noisy,
            body_dec_noisy,
            moon_ra_noisy,
            moon_dec_noisy,
        )
        min_separation_deg = float(np.min(fine_separation_deg))

        # Skip windows without partial occultation.
        if min_separation_deg > outer_threshold_deg:
            continue

        (
            external_ingress,
            internal_ingress,
            internal_egress,
            external_egress,
        ) = extract_contact_times(
            fine_days,
            fine_separation_deg,
            outer_threshold_deg,
            inner_threshold_deg,
        )

        # If complete occultation is impossible,
        # internal contacts are not meaningful.
        if not complete_occultation_possible:
            internal_ingress = None
            internal_egress = None

        external_ingress_snapshot = contact_snapshot(
            contact_points,
            fine_days,
            external_ingress,
        )
        external_egress_snapshot = contact_snapshot(
            contact_points,
            fine_days,
            external_egress,
        )

        ingress_egress_center_angle_deg = vector_angle_deg(
            center_to_body_vector(external_ingress_snapshot),
            center_to_body_vector(external_egress_snapshot),
        )

        event_rows.append(
            {
                "body": body,
                "event_index": event_index,
                "external_ingress": (
                    day_offset_to_date_string(external_ingress)
                    if external_ingress is not None
                    else ""
                ),
                "internal_ingress": (
                    day_offset_to_date_string(internal_ingress)
                    if internal_ingress is not None
                    else ""
                ),
                "internal_egress": (
                    day_offset_to_date_string(internal_egress)
                    if internal_egress is not None
                    else ""
                ),
                "external_egress": (
                    day_offset_to_date_string(external_egress)
                    if external_egress is not None
                    else ""
                ),
                "external_occultation_time": format_occultation_duration(
                    external_ingress,
                    external_egress,
                ),
                "internal_occultation_time": format_occultation_duration(
                    internal_ingress,
                    internal_egress,
                ),
                "ingress_egress_center_angle_deg": (
                    ingress_egress_center_angle_deg
                    if ingress_egress_center_angle_deg is not None
                    else ""
                ),
            }
        )
        event_index += 1

    return event_rows


def worker(body_queue, result_queue):
    moon_data = load_trajectory("moon")
    moon_apparent_radius_deg = BODY_METADATA["moon"]["apparent_radius_deg"]

    while True:
        item = body_queue.get()
        if item is None:
            body_queue.task_done()
            break

        body, _, _, body_apparent_radius_deg = item
        try:
            result = process_body(
                body,
                body_apparent_radius_deg,
                moon_data,
                moon_apparent_radius_deg,
            )
            result_queue.put(("ok", result))
        except Exception as exc:
            result_queue.put(("error", body, str(exc)))
        finally:
            body_queue.task_done()


def main():
    body_jobs = [body_job for body_job in CELESTIAL_BODIES if body_job[0] != "moon"]

    body_queue = mp.JoinableQueue()
    result_queue = mp.Queue()

    workers = []
    for _ in range(WORKER_COUNT):
        process = mp.Process(target=worker, args=(body_queue, result_queue))
        process.start()
        workers.append(process)

    for body_job in body_jobs:
        body_queue.put(body_job)

    for _ in range(WORKER_COUNT):
        body_queue.put(None)

    body_queue.join()

    all_rows = []
    errors = []
    for _ in range(len(body_jobs)):
        message = result_queue.get()
        if message[0] == "ok":
            all_rows.extend(message[1])
        else:
            _, body, err = message
            errors.append((body, err))

    for process in workers:
        process.join()

    if errors:
        for body, err in errors:
            print(f"Failed {body.upper()}: {err}")
        raise RuntimeError("One or more body jobs failed.")

    output_columns = [
        "body",
        "event_index",
        "external_ingress",
        "internal_ingress",
        "internal_egress",
        "external_egress",
        "external_occultation_time",
        "internal_occultation_time",
        "ingress_egress_center_angle_deg",
    ]

    events_df = pd.DataFrame(all_rows, columns=output_columns)
    if not events_df.empty:
        events_df = events_df.sort_values(
            by=["body", "external_ingress", "event_index"],
            kind="stable",
        ).reset_index(drop=True)

    events_df.to_csv(output_tsv_path, sep="\t", index=False)
    print(f"Wrote occultation event table to: {output_tsv_path}")


if __name__ == "__main__":
    main()
