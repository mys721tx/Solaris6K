import csv
import multiprocessing as mp
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, get_body, AltAz, solar_system_ephemeris
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.timeseries import LombScargle
from astropy.utils import iers
from scipy.signal import find_peaks

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
    "1726-01-01 00:00:00",
    scale="tt",
)

sampling_rate = 48.0  # 48 samples per day = every half hour
total_days = 365.25 * 300  # 300 Years
time_steps = np.arange(0, total_days, 1.0 / sampling_rate)
time_vector = start_time + TimeDelta(time_steps * u.day)

noise_std = 2.0 / 60
np.random.seed(42)

# Keep body, x-limit, and known orbital period (days) together so entries can be toggled in one place.
celestial_bodies = [
    ("sun", 0.007, 365.25),  # Sidereal year
    ("moon", 0.08, 27.322),  # Sidereal month
    ("mercury", 0.02, 87.969),
    ("venus", 0.005, 224.701),
    ("mars", 0.003, 686.971),
    ("jupiter", 0.003, 4332.589),
    ("saturn", 0.003, 10759.22),
]

# Dense frequency grid focused on the high-value astronomical zones
min_freq = 1.0 / total_days
max_freq = 2
frequencies = np.linspace(min_freq, max_freq, 1000000)

print("Pre-calculating observation framework...")
altaz_frame = AltAz(obstime=time_vector, location=location)

fft_output_dir = Path("fft_tables")
plot_output_dir = Path("fft_plots")


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
    # This allows the algorithm to immediately drop into a valley and find a second peak.
    safety_distance_bins = 15

    # 2. Advanced Peak Finding relying purely on local topographic prominence
    peaks, _ = find_peaks(
        combined_power,
        prominence=0.03
        * np.max(combined_power),  # Captures distinct features > 3% of max power
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

    with fft_table_path.open(mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "cycle_per_day",
                "period_days",
                "period_hours",
                "power_ra",
                "power_dec",
                "power_combined",
            ]
        )
        for i in range(len(top_peaks)):
            writer.writerow(
                [
                    f"{top_freqs[i]:.5f}",
                    f"{top_period_days[i]:.5f}",
                    f"{top_period_hours[i]:.5f}",
                    f"{top_power_ra[i]:.5f}",
                    f"{top_power_dec[i]:.5f}",
                    f"{top_power_combined[i]:.5f}",
                ]
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
    start_date_label = str(start_time).split()[0]
    ax1.scatter(
        ra_detrended[first_orbit_mask],
        dec_detrended[first_orbit_mask],
        s=2,
        color="crimson",
        alpha=0.5,
        label=f"First orbit from {start_date_label}",
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

        body, x_limit, orbital_period = item
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

    worker_count = 8
    body_queue = mp.JoinableQueue()
    result_queue = mp.Queue()

    workers = []
    for _ in range(worker_count):
        process = mp.Process(target=worker, args=(body_queue, result_queue))
        process.start()
        workers.append(process)

    for body_job in celestial_bodies:
        body_queue.put(body_job)

    for _ in range(worker_count):
        body_queue.put(None)

    body_queue.join()

    errors = []
    results = []
    for _ in range(len(celestial_bodies)):
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
