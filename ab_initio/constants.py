TOTAL_DAYS = 365.25 * 300  # 300 Years

START_TIME_TT = "1726-01-01 00:00:00"
START_DATE_LABEL = START_TIME_TT.split()[0]

# Keep body, x-limit, and known orbital period (days) together
# so entries can be toggled in one place.
CELESTIAL_BODIES = [
    ("sun", 0.007, 365.25),  # Sidereal year
    ("moon", 0.08, 27.322),  # Sidereal month
    ("mercury", 0.02, 87.969),
    ("venus", 0.005, 224.701),
    ("mars", 0.003, 686.971),
    ("jupiter", 0.003, 4332.589),
    ("saturn", 0.003, 10759.22),
]

TRAJECTORY_TABLES_DIRNAME = "trajectories"

WORKER_COUNT = 8
