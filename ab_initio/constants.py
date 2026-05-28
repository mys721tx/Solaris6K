TOTAL_DAYS = 365.25 * 300  # 300 Years

START_TIME_TT = "1726-01-01 00:00:00"
START_DATE_LABEL = START_TIME_TT.split()[0]

# Keep body, x-limit, known orbital period (days), and representative
# apparent angular radius (degrees) together
# so entries can be toggled in one place.
CELESTIAL_BODIES = [
    ("sun", 0.007, 365.25, 0.2666),  # Sidereal year
    ("moon", 0.08, 27.322, 0.2725),  # Sidereal month
    ("mercury", 0.02, 87.969, 0.0012),
    ("venus", 0.005, 224.701, 0.0040),
    ("mars", 0.003, 686.971, 0.0018),
    ("jupiter", 0.003, 4332.589, 0.0055),
    ("saturn", 0.003, 10759.22, 0.0024),
]

TRAJECTORY_TABLES_DIRNAME = "trajectories"

WORKER_COUNT = 8
