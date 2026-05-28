from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import (
    CELESTIAL_BODIES,
    START_TIME_TT,
    TRAJECTORY_TABLES_DIRNAME,
)

input_dir = Path(TRAJECTORY_TABLES_DIRNAME)
output_plot_path = Path("moon_occultation_instances.png")
SIMULATION_START = datetime.strptime(START_TIME_TT, "%Y-%m-%d %H:%M:%S")

BODY_METADATA = {
    body: {
        "x_limit": x_limit,
        "orbital_period": orbital_period,
        "apparent_radius_deg": apparent_radius_deg,
    }
    for body, x_limit, orbital_period, apparent_radius_deg in CELESTIAL_BODIES
}


def load_trajectory(body):
    trajectory_hdf_path = input_dir / f"{body}_trajectory.h5"
    if not trajectory_hdf_path.exists():
        raise FileNotFoundError(f"Missing input HDF5 for {body}: {trajectory_hdf_path}")

    trajectory_data = pd.read_hdf(
        trajectory_hdf_path,
        key="trajectory",
    )

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


def detect_occultation_instances(
    moon_data,
    body_data,
    moon_apparent_radius_deg,
    body_apparent_radius_deg,
):
    merged = moon_data.merge(
        body_data,
        on="time_point",
        how="inner",
        suffixes=("_moon", "_body"),
    )

    if merged.empty:
        return merged, merged

    merged["separation_deg"] = angular_separation_deg(
        merged["ra_filtered_moon"],
        merged["dec_filtered_moon"],
        merged["ra_filtered_body"],
        merged["dec_filtered_body"],
    )

    # Complete occultation: the full angular disk of the occulted body
    # is contained inside the Moon's apparent disk.
    merged["is_complete_occultation"] = (
        merged["separation_deg"] + body_apparent_radius_deg
    ) <= moon_apparent_radius_deg

    occultations = merged.loc[
        merged["is_complete_occultation"],
        ["time_point", "separation_deg"],
    ].copy()
    return merged, occultations


def extract_complete_occultation_events(merged_data):
    if merged_data.empty:
        return []

    times = merged_data["time_point"].to_numpy(dtype=float)
    is_complete = merged_data["is_complete_occultation"].to_numpy(dtype=bool)
    if not np.any(is_complete):
        return []

    time_gaps = np.diff(times)
    valid_gaps = time_gaps[time_gaps > 0.0]
    if len(valid_gaps) == 0:
        nominal_dt = 0.0
    else:
        nominal_dt = float(np.median(valid_gaps))
    split_threshold = nominal_dt * 1.5

    events = []
    in_event = False
    start_index = 0

    for index in range(len(times)):
        if not in_event and is_complete[index]:
            in_event = True
            start_index = index

        if not in_event:
            continue

        is_last_sample = index == len(times) - 1
        next_breaks_event = False
        if not is_last_sample:
            next_breaks_event = (not is_complete[index + 1]) or (
                times[index + 1] - times[index] > split_threshold
            )

        if is_last_sample or next_breaks_event:
            end_index = index
            start_time = float(times[start_index])
            end_time = float(times[end_index])
            sample_count = int(end_index - start_index + 1)
            duration_days = max(0.0, end_time - start_time)
            events.append(
                {
                    "start_time": start_time,
                    "end_time": end_time,
                    "duration_days": duration_days,
                    "sample_count": sample_count,
                }
            )
            in_event = False

    return events


def plot_occultation_timeline(occultation_events, plot_path):
    if not occultation_events:
        print("No Moon occultation instances found for any body.")
        return

    bodies = sorted(occultation_events.keys())
    body_to_y = {body: index for index, body in enumerate(bodies)}

    fig_height = max(4.5, 0.6 * len(bodies) + 2.0)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    for body in bodies:
        events = occultation_events[body]
        ax.scatter(
            events["time_point"],
            np.full(len(events), body_to_y[body]),
            s=12,
            alpha=0.8,
            label=body.upper(),
        )

    ax.set_yticks(list(body_to_y.values()))
    ax.set_yticklabels([body.upper() for body in bodies])
    ax.set_xlabel("Simulation Time (days from start)")
    ax.set_ylabel("Occulted Body")
    ax.set_title("Moon Occultation Instances")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", ncol=2, fontsize=8)

    fig.tight_layout()
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    print(f"Saved occultation timeline plot: {plot_path}")


def day_offset_to_date_string(day_offset):
    timestamp = SIMULATION_START + timedelta(days=float(day_offset))
    return timestamp.strftime("%Y-%m-%d %H:%M:%S")


def summarize_occultations(body, merged_data, occultations, events):
    sample_count = len(merged_data)
    occultation_count = len(occultations)
    if occultation_count == 0:
        print(
            f"{body.upper()}: 0 complete occultation samples out of {sample_count} "
            "shared samples"
        )
        return

    first_time = float(occultations["time_point"].iloc[0])
    last_time = float(occultations["time_point"].iloc[-1])
    min_sep = float(occultations["separation_deg"].min())
    first_date = day_offset_to_date_string(first_time)
    last_date = day_offset_to_date_string(last_time)

    print(
        f"{body.upper()}: {occultation_count} complete occultation samples out of "
        f"{sample_count} shared samples; first={first_date}, "
        f"last={last_date}, "
        f"min separation={min_sep:.4f} deg"
    )

    total_duration_days = sum(event["duration_days"] for event in events)
    print(
        f"{body.upper()}: {len(events)} complete occultation events, "
        f"total complete-occultation duration={total_duration_days:.4f} days"
    )
    for index, event in enumerate(events, start=1):
        start_date = day_offset_to_date_string(event["start_time"])
        end_date = day_offset_to_date_string(event["end_time"])
        print(
            f"  Event {index}: start={start_date}, "
            f"end={end_date}, "
            f"duration={event['duration_days']:.6f} days, "
            f"samples={event['sample_count']}"
        )


def main():
    moon_apparent_radius_deg = BODY_METADATA["moon"]["apparent_radius_deg"]
    moon_data = load_trajectory("moon")
    occultation_events = {}

    for body, _, _, body_apparent_radius_deg in CELESTIAL_BODIES:
        if body == "moon":
            continue

        body_data = load_trajectory(body)
        merged_data, occultations = detect_occultation_instances(
            moon_data,
            body_data,
            moon_apparent_radius_deg,
            body_apparent_radius_deg,
        )
        events = extract_complete_occultation_events(merged_data)
        summarize_occultations(body, merged_data, occultations, events)

        if not occultations.empty:
            occultation_events[body] = occultations

    plot_occultation_timeline(occultation_events, output_plot_path)


if __name__ == "__main__":
    main()
