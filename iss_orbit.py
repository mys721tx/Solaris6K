# %%
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos
from datetime import datetime, timezone, timedelta

# %%
# --- Parameters ---
tle_url = "https://celestrak.org/NORAD/elements/stations.txt"
observer_lat, observer_lon, observer_elev = (
    34.41916310,
    -119.87076250,
    20,
)  #
duration_hours = 24
dt_seconds = 60
altitude_cut_deg = 10  # Horizon cutoff - set values to NaN below this

# --- Load ISS TLE and Observer ---
print("Downloading TLE...")
satellites = load.tle_file(tle_url)
iss = {sat.name: sat for sat in satellites}["ISS (ZARYA)"]
ts = load.timescale()
observer = Topos(
    latitude_degrees=observer_lat,
    longitude_degrees=observer_lon,
    elevation_m=observer_elev,
)

# --- Generate time samples ---
total_steps = int(duration_hours * 3600 / dt_seconds)
now_utc = datetime.now(timezone.utc)
times = ts.utc(
    now_utc.year,
    now_utc.month,
    now_utc.day,
    hour=now_utc.hour,
    minute=now_utc.minute + np.arange(total_steps) * dt_seconds / 60.0,
)

# --- Calculate azimuth and elevation ---
difference = iss - observer
topo = difference.at(times)
alt, az, _ = topo.altaz()
alt_deg = alt.degrees
az_deg_raw = az.degrees

# --- Prepare data for plotting and speed calculations ---
elev_plot = np.array(alt_deg, dtype=float)  # Keep all values for plotting
az_plot = np.unwrap(np.radians(np.array(az_deg_raw)))  # Keep all values for plotting

# For speed calculations, set below-horizon values to NaN
elev_tracked = np.array(alt_deg, dtype=float)
az_tracked = np.unwrap(np.radians(np.array(az_deg_raw)))
elev_tracked[elev_tracked <= altitude_cut_deg] = np.nan
az_tracked[np.array(alt_deg) <= altitude_cut_deg] = np.nan

# --- Convert time to seconds and hours for plotting ---
elapsed_sec = (times - times[0]) * 86400.0
time_hours = elapsed_sec / 3600.0

# --- Compute angular velocities (deg/sec) ---
d_elev = np.gradient(elev_tracked, elapsed_sec)
d_az = np.gradient(az_tracked, elapsed_sec)  # rad/sec
d_az_deg = np.degrees(d_az)  # deg/sec

# --- Detect flybys (periods above cutoff elevation) ---
above_cutoff = elev_plot > altitude_cut_deg
# Find transitions
transitions = np.diff(above_cutoff.astype(int))
flyby_starts = np.where(transitions == 1)[0] + 1  # Rising edge
flyby_ends = np.where(transitions == -1)[0] + 1  # Falling edge

# Handle edge cases
if above_cutoff[0]:
    flyby_starts = np.concatenate([[0], flyby_starts])
if above_cutoff[-1]:
    flyby_ends = np.concatenate([flyby_ends, [len(above_cutoff) - 1]])

# Calculate flyby durations and mid-points
flyby_durations = []
flyby_midtimes = []
for start, end in zip(flyby_starts, flyby_ends):
    duration_sec = elapsed_sec[end] - elapsed_sec[start]
    duration_min = duration_sec / 60.0
    midtime = (time_hours[start] + time_hours[end]) / 2.0
    flyby_durations.append(duration_min)
    flyby_midtimes.append(midtime)

# --- Plotting ---
fig, axs = plt.subplots(
    3, 1, figsize=(12, 13), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1]}
)

# Top panel: Elevation
axs[0].plot(time_hours, elev_plot, "b-", linewidth=1.5)
axs[0].axhline(
    y=altitude_cut_deg,
    color="r",
    linestyle="--",
    linewidth=1,
    alpha=0.7,
    label=f"Cutoff ({altitude_cut_deg}°)",
)
axs[0].set_ylabel("Elevation [deg]")
axs[0].set_title(f"ISS Tracking at {observer_lat:.2f}°, {observer_lon:.2f}°")
axs[0].legend()
axs[0].grid(True, alpha=0.3)

# Middle panel: Azimuth
axs[1].plot(time_hours, np.degrees(az_plot) % 360, "r-", linewidth=1.5)
axs[1].set_ylabel("Azimuth [deg]")
axs[1].grid(True, alpha=0.3)

# Bottom panel: Angular velocities
axs[2].plot(time_hours, d_elev, label="d(El)/dt", linewidth=1.5)
axs[2].plot(time_hours, d_az_deg, label="d(Az)/dt", linewidth=1.5)
axs[2].set_ylabel("Angular Velocity [deg/s]")
axs[2].set_xlabel("Time [hours]")
axs[2].legend()
axs[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %%
