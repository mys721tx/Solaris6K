import multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import (
    AltAz,
    EarthLocation,
    get_body,
    solar_system_ephemeris,
)
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.utils import iers

from constants import (
    CELESTIAL_BODIES,
    START_TIME_TT,
    TOTAL_DAYS,
    TRAJECTORY_TABLES_DIRNAME,
    WORKER_COUNT,
)

# Tell astropy to silently ignore the missing historical EOP data
iers.conf.iers_degraded_accuracy = "ignore"

# ==========================================
# 1. Setup Simulation Constants
# ==========================================

# Use the local DE440 ephemeris (covers 1550–2650 AD).
# Switch to de421.bsp for broader historical range (−3000 to +3000 AD).
solar_system_ephemeris.set("de440")

location = EarthLocation(
    lat=0 * u.deg,
    lon=0 * u.deg,
    height=0 * u.m,
)
# Use scale="tt" (Terrestrial Time) for pre-1960 dates.
# UTC is only formally defined from 1960 onward.
start_time = Time(
    START_TIME_TT,
    scale="tt",
)

sampling_rate = 48.0  # 48 samples per day = every 30 minutes
time_steps = np.arange(0, TOTAL_DAYS, 1.0 / sampling_rate)
time_vector = start_time + TimeDelta(time_steps * u.day)

noise_std = 2.0 / 60
np.random.seed(42)

print("Pre-calculating observation framework...")
altaz_frame = AltAz(obstime=time_vector, location=location)

output_dir = Path(TRAJECTORY_TABLES_DIRNAME)


def process_body(body, x_limit, orbital_period):
    print(f"Analyzing {body.upper()}...")

    # Use deterministic body-specific noise so multiprocessing
    # output is stable.
    rng = np.random.default_rng(42 + sum(ord(c) for c in body))

    # Fetch and inject noise
    coords = get_body(
        body,
        time_vector,
        location=location,
    )
    altaz = coords.transform_to(altaz_frame)
    alt_pure = altaz.alt.deg

    ra_pure = coords.ra.deg
    dec_pure = coords.dec.deg

    ra_noisy = ra_pure + rng.normal(
        0,
        noise_std,
        size=len(ra_pure),
    )
    dec_noisy = dec_pure + rng.normal(
        0,
        noise_std,
        size=len(dec_pure),
    )

    # Filter by horizon
    horizon_mask = alt_pure > 0
    filtered_days = time_steps[horizon_mask]
    ra_filtered = ra_noisy[horizon_mask]
    dec_filtered = dec_noisy[horizon_mask]

    # Detrend
    ra_unwrapped = np.unwrap(np.radians(ra_filtered))
    dec_radians = np.radians(dec_filtered)

    poly_ra = np.polyfit(filtered_days, ra_unwrapped, 1)
    ra_detrended = np.degrees(
        ra_unwrapped - np.polyval(poly_ra, filtered_days),
    )

    poly_dec = np.polyfit(filtered_days, dec_radians, 1)
    dec_detrended = np.degrees(
        dec_radians
        - np.polyval(
            poly_dec,
            filtered_days,
        )
    )

    trajectory_hdf_path = output_dir / f"{body}_trajectory.h5"
    pd.DataFrame(
        {
            "time_point": filtered_days,
            "ra_filtered": ra_filtered,
            "dec_filtered": dec_filtered,
            "ra_detrended": ra_detrended,
            "dec_detrended": dec_detrended,
        }
    ).to_hdf(
        trajectory_hdf_path,
        key="trajectory",
        mode="w",
    )

    return {
        "body": body,
        "trajectory_hdf_path": trajectory_hdf_path,
    }


def worker(body_queue, result_queue):
    while True:
        item = body_queue.get()
        if item is None:
            body_queue.task_done()
            break

        body, x_limit, orbital_period, _ = item
        try:
            result = process_body(body, x_limit, orbital_period)
            result_queue.put(("ok", result))
        except Exception as exc:
            result_queue.put(("error", body, str(exc)))
        finally:
            body_queue.task_done()


def main():
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing windowed coordinate tables to: {output_dir}")

    print("\n=== STARTING BULK SIGNAL PROCESSING COHORT ===")

    body_queue = mp.JoinableQueue()
    result_queue = mp.Queue()

    workers = []
    for _ in range(WORKER_COUNT):
        process = mp.Process(target=worker, args=(body_queue, result_queue))
        process.start()
        workers.append(process)

    for body_job in CELESTIAL_BODIES:
        body_queue.put(body_job)

    for _ in range(WORKER_COUNT):
        body_queue.put(None)

    body_queue.join()

    errors = []
    results = []
    for _ in range(len(CELESTIAL_BODIES)):
        message = result_queue.get()
        if message[0] == "ok":
            results.append(message[1])
        else:
            _, body, err = message
            errors.append((body, err))

    for process in workers:
        process.join()

    for result in sorted(results, key=lambda item: item["body"]):
        body = result["body"]
        print(
            f"Saved {body.upper()} trajectory data: " f"{result['trajectory_hdf_path']}"
        )

    if errors:
        for body, err in errors:
            print(f"Failed {body.upper()}: {err}")
        raise RuntimeError("One or more body jobs failed.")

    print(f"All trajectory data written to: {output_dir}")


if __name__ == "__main__":
    main()
