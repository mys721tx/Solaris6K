import multiprocessing as mp
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from constants import (
    CELESTIAL_BODIES,
    START_DATE_LABEL,
    TOTAL_DAYS,
    TRAJECTORY_TABLES_DIRNAME,
    WORKER_COUNT,
)

# Dense frequency grid focused on the high-value astronomical zones
min_freq = 1.0 / TOTAL_DAYS
max_freq = 2
frequencies = np.linspace(min_freq, max_freq, 1000000)

print("Pre-calculating observation framework...")

fft_output_dir = Path("fft_tables")
plot_output_dir = Path("fft_plots")
input_dir = Path(TRAJECTORY_TABLES_DIRNAME)


def process_body(body, x_limit, orbital_period):
    print(f"Analyzing {body.upper()}...")

    trajectory_hdf_path = input_dir / f"{body}_trajectory.h5"
    if not trajectory_hdf_path.exists():
        raise FileNotFoundError(f"Missing input HDF5 for {body}: {trajectory_hdf_path}")

    trajectory_data = pd.read_hdf(
        trajectory_hdf_path,
        key="trajectory",
    )

    filtered_days = trajectory_data["time_point"]
    ra_detrended = trajectory_data["ra_detrended"]
    dec_detrended = trajectory_data["dec_detrended"]

    if trajectory_data.empty:
        raise ValueError(f"Input HDF5 for {body} is empty: {trajectory_hdf_path}")

    # Window detrended signal prior to spectral estimation to reduce leakage.
    window = np.hanning(len(filtered_days))
    ra_windowed = ra_detrended * window
    dec_windowed = dec_detrended * window

    # Compute Lomb-Scargle
    power_ra = LombScargle(
        filtered_days,
        ra_windowed,
    ).power(
        frequencies,
    )
    power_dec = LombScargle(
        filtered_days,
        dec_windowed,
    ).power(
        frequencies,
    )
    combined_power = power_ra + power_dec

    fft_table_path = fft_output_dir / f"{body}_fft.tsv"

    # 1. Very aggressive distance gate (e.g., just 10-20 frequency bins)
    # This allows the algorithm to immediately drop into a valley and
    # find a second peak.
    safety_distance_bins = 15

    # 2. Advanced Peak Finding relying purely on local topographic prominence
    peaks, _ = find_peaks(
        combined_power,
        prominence=0.03 * np.max(combined_power),
        distance=safety_distance_bins,
    )

    # Write top 20 peaks (by combined power) to TSV.
    top_peaks = peaks[np.argsort(combined_power[peaks])[::-1]][:20]
    top_freqs = frequencies[top_peaks]
    top_period_days = 1.0 / top_freqs
    top_period_hours = top_period_days * 24.0
    top_power_ra = power_ra[top_peaks]
    top_power_dec = power_dec[top_peaks]
    top_power_combined = combined_power[top_peaks]

    pd.DataFrame(
        {
            "cycle_per_day": top_freqs,
            "period_days": top_period_days,
            "period_hours": top_period_hours,
            "power_ra": top_power_ra,
            "power_dec": top_power_dec,
            "power_combined": top_power_combined,
        }
    ).to_csv(
        fft_table_path,
        sep="\t",
        index=False,
        float_format="%.5f",
    )

    visible_peaks = peaks[(frequencies[peaks] >= 0.0) & (frequencies[peaks] <= x_limit)]
    if len(visible_peaks) > 0:
        sorted_peaks = visible_peaks[np.argsort(combined_power[visible_peaks])[::-1]]
    else:
        sorted_peaks = np.array([], dtype=int)

    # Trajectory
    fig1, ax1 = plt.subplots(figsize=(8, 6))
    ax1.scatter(
        ra_detrended,
        dec_detrended,
        s=0.5,
        color="gray",
        alpha=0.05,
    )
    first_orbit_mask = filtered_days < orbital_period
    ax1.scatter(
        ra_detrended[first_orbit_mask],
        dec_detrended[first_orbit_mask],
        s=2,
        color="crimson",
        alpha=0.5,
        label=f"First orbit from {START_DATE_LABEL}",
    )
    ax1.set_title(f"2D Celestial Trajectory: {body.upper()}")
    ax1.set_xlabel("Relative Right Ascension (Degrees)")
    ax1.set_ylabel("Relative Declination (Degrees)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right")
    fig1.tight_layout()
    trajectory_path = plot_output_dir / f"{body}_trajectory.png"
    fig1.savefig(trajectory_path, dpi=300)
    plt.close(fig1)

    # Power spectrum
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    ax2.plot(frequencies, combined_power, color="midnightblue", lw=1.2)

    # Apply the per-body viewport from the tuple config.
    ax2.set_xlim(0.0, x_limit)

    # Annotate the top 5 dominant landmarks within the visual limits
    for i in range(min(5, len(sorted_peaks))):
        p_idx = int(sorted_peaks[i])
        p_freq = float(frequencies[p_idx])
        p_days = 1.0 / p_freq
        p_label = f"{p_days:.2f}d" if p_days >= 1.0 else f"{p_days*24:.1f}h"
        ax2.axvline(
            p_freq,
            color="crimson" if i == 0 else "forestgreen",
            linestyle=":" if i == 0 else "--",
            alpha=0.6,
        )
        ax2.text(
            p_freq,
            max(combined_power) * (0.85 - i * 0.08),
            f" {p_label}",
            color="crimson" if i == 0 else "forestgreen",
            fontsize=9,
            weight="bold",
        )

    ax2.set_title(f"Lomb-Scargle Spectrum: {body.upper()}")
    ax2.set_xlabel("Cycles / Day")
    ax2.set_ylabel("Spectral Power")
    ax2.grid(True, linestyle="--", alpha=0.5)
    fig2.tight_layout()
    spectrum_path = plot_output_dir / f"{body}_spectrum.png"
    fig2.savefig(spectrum_path, dpi=300)
    plt.close(fig2)

    return {
        "body": body,
        "fft_table_path": fft_table_path,
        "trajectory_path": trajectory_path,
        "spectrum_path": spectrum_path,
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
    fft_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing high-precision FFT tables to: {fft_output_dir}")

    plot_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing plot outputs to: {plot_output_dir}")

    print("\n=== STARTING BULK SIGNAL PROCESSING + PLOTTING COHORT ===")

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
        print(f"Saved FFT table: {result['fft_table_path']}")
        print(
            f"Saved {body.upper()} plots: "
            f"{result['trajectory_path']}, {result['spectrum_path']}"
        )

    if errors:
        for body, err in errors:
            print(f"Failed {body.upper()}: {err}")
        raise RuntimeError("One or more body jobs failed.")

    print(f"All subplot files written to: {plot_output_dir}")


if __name__ == "__main__":
    main()
