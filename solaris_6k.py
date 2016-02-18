"""
Solaris6K.py: Calculate approximate positions of the major planets
"""

import numpy as np

import constants


def sin_wrapper(theta):
    """
    Warpper for np.sin to take degree
    """

    return np.sin(np.deg2rad(theta))


def arcsin_wrapper(num):
    """
    Warpper for np.arcsin to return degree
    """

    return np.rad2deg(np.arcsin(num))


def cos_wrapper(theta):
    """
    Warpper for np.cos to take degree
    """

    return np.cos(np.deg2rad(theta))


def arccos_wrapper(num):
    """
    Warpper for np.arccos to return degree
    """

    return np.rad2deg(np.arccos(num))


def get_century(jd_time):
    """
    Take a JD time and return the century pased since J2000.0
    """

    num_century = (jd_time - constants.J2000_0) / constants.CONVERSION_FACTOR

    if abs(num_century) > 30:
        raise ValueError("Date is outside of scope (3000 BCE to 3000 CE).")
    else:
        return num_century


def get_current_value(base_value, rate, time):
    """
    Calcuate the value at given base_value, rate, and time.
    """

    return base_value + rate * time


def get_perihelion_argument(
    perihelion_longitude,
    ascending_node_longitude,
):
    """
    Calculate the argument of perihelion.
    """

    return perihelion_longitude - ascending_node_longitude


def get_mean_anomaly(
    mean_longitude,
    perihelion_longitude,
    additional_terms,
    time,
):
    """
    Calculate the mean anomaly and normalize it within +/- 180 degree
    """

    term_b, term_c, term_s, term_f = additional_terms

    mean_anomaly = mean_longitude - perihelion_longitude

    # Additional terms for planets after Jupiter
    mean_anomaly += term_b * time**2
    mean_anomaly += term_c * np.cos(term_f * time)
    mean_anomaly += term_s * np.sin(term_f * time)

    # Normalize the mean_anomaly
    return arcsin_wrapper(sin_wrapper(mean_anomaly))


def get_eccentric_anomaly(
    mean_anomaly,
    eccentricity,
    tol,
):
    """
    Iterate solver for eccentric anomaly
    """

    term_e = eccentricity * constants.ECCENTRICITIY_FACTOR
    eccentric_anomaly = mean_anomaly + term_e * sin_wrapper(mean_anomaly)
    eccentric_anomaly_delta = np.float64("inf")

    while any(eccentric_anomaly_delta > tol):
        mean_anomaly_delta = mean_anomaly - (
            eccentric_anomaly - term_e * sin_wrapper(eccentric_anomaly)
        )
        eccentric_anomaly_delta = mean_anomaly_delta / (
            1 - eccentricity * cos_wrapper(eccentric_anomaly)
        )
        eccentric_anomaly += eccentric_anomaly_delta
    return eccentric_anomaly


def get_ecliptic_coordinate(
    heliocentric_coordinate,
    perihelion_argument,
    ascending_node_longitude,
    inclination,
):
    """
    Convert heliocentric coordinate to J2000 ecliptic plane coordinate
    """

    ecliptic_x_matrix = np.column_stack(
        (
            cos_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude)
            - sin_wrapper(perihelion_argument)
            * sin_wrapper(ascending_node_longitude)
            * cos_wrapper(inclination),
            -sin_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude)
            - cos_wrapper(perihelion_argument)
            * sin_wrapper(ascending_node_longitude)
            * cos_wrapper(inclination),
            np.array(
                [0 for num in enumerate(perihelion_argument)],
                np.dtype("float64"),
            ),
        )
    )

    ecliptic_y_matrix = np.column_stack(
        (
            cos_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude)
            + sin_wrapper(perihelion_argument)
            * cos_wrapper(ascending_node_longitude)
            * cos_wrapper(inclination),
            -sin_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude)
            + cos_wrapper(perihelion_argument)
            * cos_wrapper(ascending_node_longitude)
            * cos_wrapper(inclination),
            np.array(
                [0 for num in enumerate(perihelion_argument)],
                np.dtype("float64"),
            ),
        )
    )

    ecliptic_z_matrix = np.column_stack(
        (
            sin_wrapper(perihelion_argument) * sin_wrapper(inclination),
            cos_wrapper(perihelion_argument) * sin_wrapper(inclination),
            np.array(
                [0 for num in enumerate(perihelion_argument)],
                np.dtype("float64"),
            ),
        )
    )

    return np.column_stack(
        (
            np.sum(heliocentric_coordinate * ecliptic_x_matrix, axis=1),
            np.sum(heliocentric_coordinate * ecliptic_y_matrix, axis=1),
            np.sum(heliocentric_coordinate * ecliptic_z_matrix, axis=1),
        )
    )


def main():
    """
    Calculations
    """

    time = 2451545.0
    century_since_epoch = get_century(time)

    semi_major_axes = get_current_value(
        constants.SEMI_MAJOR_AXES,
        constants.SEMI_MAJOR_AXES_RATE,
        century_since_epoch,
    )

    eccentricities = get_current_value(
        constants.ECCENTRICITIES,
        constants.ECCENTRICITIES_RATE,
        century_since_epoch,
    )

    inclinations = get_current_value(
        constants.INCLINATIONS,
        constants.INCLINATIONS_RATE,
        century_since_epoch,
    )

    mean_longitudes = get_current_value(
        constants.MEAN_LONGITUDES,
        constants.MEAN_LONGITUDES_RATE,
        century_since_epoch,
    )

    perihelion_longitudes = get_current_value(
        constants.PERIHELION_LONGITUDES,
        constants.PERIHELION_LONGITUDES_RATE,
        century_since_epoch,
    )

    ascending_node_longitudes = get_current_value(
        constants.ASCENDING_NODE_LONGITUDES,
        constants.ASCENDING_NODE_LONGITUDES_RATE,
        century_since_epoch,
    )

    perihelion_arguments = get_perihelion_argument(
        perihelion_longitudes, ascending_node_longitudes
    )

    mean_anomalies = get_mean_anomaly(
        mean_longitudes,
        perihelion_longitudes,
        (
            constants.TERM_B,
            constants.TERM_C,
            constants.TERM_S,
            constants.TERM_F,
        ),
        century_since_epoch,
    )

    eccentric_anomalies = get_eccentric_anomaly(
        mean_anomalies,
        eccentricities,
        np.array(
            [constants.TOL for num in enumerate(mean_anomalies)],
        ),
    )

    heliocentric_coordinates = np.column_stack(
        (
            semi_major_axes * (cos_wrapper(eccentric_anomalies) - eccentricities),
            semi_major_axes
            * np.sqrt(1 - eccentricities**2)
            * sin_wrapper(eccentric_anomalies),
            np.array(
                [0 for num in enumerate(mean_anomalies)],
                np.dtype("float64"),
            ),
        )
    )

    ecliptic_coordinates = get_ecliptic_coordinate(
        heliocentric_coordinates,
        perihelion_arguments,
        ascending_node_longitudes,
        inclinations,
    )

    eq_matrix = np.matrix(
        [
            [1, 0, 0],
            [
                0,
                cos_wrapper(constants.OBLIQUITY),
                -sin_wrapper(constants.OBLIQUITY),
            ],
            [
                0,
                sin_wrapper(constants.OBLIQUITY),
                cos_wrapper(constants.OBLIQUITY),
            ],
        ]
    )

    print(np.dot(eq_matrix, ecliptic_coordinates.T).T)


if __name__ == "__main__":
    main()
