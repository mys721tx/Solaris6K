"""
Solaris6K.py: Calculate approximate positions of the major planets
"""

import numpy
import Constants

def sin_wrapper(theta):
    """
    Warpper for numpy.sin to take degree
    """

    return numpy.sin(numpy.deg2rad(theta))

def arcsin_wrapper(num):
    """
    Warpper for numpy.arcsin to return degree
    """

    return numpy.rad2deg(numpy.arcsin(num))

def cos_wrapper(theta):
    """
    Warpper for numpy.cos to take degree
    """

    return numpy.cos(numpy.deg2rad(theta))

def arccos_wrapper(num):
    """
    Warpper for numpy.arccos to return degree
    """

    return numpy.rad2deg(numpy.arccos(num))

def get_century(jd_time):
    """
    Take a JD time and return the century pased since J2000.0
    """

    num_century = (jd_time - Constants.J2000_0) / Constants.CONVERSION_FACTOR

    if abs(num_century) > 30:
        raise ValueError("Date is outside of scope (3000 BCE to 3000 CE).")
    else:
        return num_century

def get_current_value(base_value, rate, time):
    """
    Calcuate the value at given base_value, rate, and time.
    """

    return base_value + rate * time

def get_perihelion_argument(perihelion_longitude, ascending_node_longitude):
    """
    Calculate the argument of perihelion.
    """

    return perihelion_longitude - ascending_node_longitude

def get_mean_anomaly(
        mean_longitude,
        perihelion_longitude,
        additional_terms,
        time
    ):
    """
    Calculate the mean anomaly and normalize it within +/- 180 degree
    """

    term_b, term_c, term_s, term_f = additional_terms

    mean_anomaly = mean_longitude - perihelion_longitude

    # Additional terms for planets after Jupiter
    mean_anomaly += term_b * time ** 2
    mean_anomaly += term_c * numpy.cos(term_f * time)
    mean_anomaly += term_s * numpy.sin(term_f * time)

    # Normalize the mean_anomaly
    return arcsin_wrapper(sin_wrapper(mean_anomaly))

def get_eccentric_anomaly(mean_anomaly, eccentricity, tol):
    """
    Iterate solver for eccentric anomaly
    """

    term_e = eccentricity * Constants.ECCENTRICITIY_FACTOR
    eccentric_anomaly = mean_anomaly + term_e * sin_wrapper(mean_anomaly)
    eccentric_anomaly_delta = numpy.float64("inf")

    while any(eccentric_anomaly_delta > tol):
        mean_anomaly_delta = mean_anomaly - (eccentric_anomaly - term_e * sin_wrapper(eccentric_anomaly))
        eccentric_anomaly_delta = mean_anomaly_delta / (1 - eccentricity * cos_wrapper(eccentric_anomaly))
        eccentric_anomaly += eccentric_anomaly_delta
    return eccentric_anomaly

def get_ecliptic_coordinate(
        heliocentric_coordinate,
        perihelion_argument,
        ascending_node_longitude,
        inclination
    ):
    """
    Convert heliocentric coordinate to J2000 ecliptic plane coordinate
    """

    ecliptic_x_matrix = numpy.column_stack(
        (
            cos_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude) - sin_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude) * cos_wrapper(inclination),
            - sin_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude) - cos_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude) * cos_wrapper(inclination),
            numpy.array(
                [0 for num in enumerate(perihelion_argument)],
                numpy.dtype("float64")
            )
        )
    )

    ecliptic_y_matrix = numpy.column_stack(
        (
            cos_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude) + sin_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude) * cos_wrapper(inclination),
            - sin_wrapper(perihelion_argument) * sin_wrapper(ascending_node_longitude) + cos_wrapper(perihelion_argument) * cos_wrapper(ascending_node_longitude) * cos_wrapper(inclination),
            numpy.array(
                [0 for num in enumerate(perihelion_argument)],
                numpy.dtype("float64")
            )
        )
    )

    ecliptic_z_matrix = numpy.column_stack(
        (
            sin_wrapper(perihelion_argument) * sin_wrapper(inclination),
            cos_wrapper(perihelion_argument) * sin_wrapper(inclination),
            numpy.array(
                [0 for num in enumerate(perihelion_argument)],
                numpy.dtype("float64")
            )
        )
    )

    return numpy.column_stack(
        (
            numpy.sum(heliocentric_coordinate * ecliptic_x_matrix, axis=1),
            numpy.sum(heliocentric_coordinate * ecliptic_y_matrix, axis=1),
            numpy.sum(heliocentric_coordinate * ecliptic_z_matrix, axis=1)
        )
    )

def main():
    """
    Calculations
    """

    time = 2451545.0
    century_since_epoch = get_century(time)

    semi_major_axes = get_current_value(
        Constants.SEMI_MAJOR_AXES,
        Constants.SEMI_MAJOR_AXES_RATE,
        century_since_epoch
    )

    eccentricities = get_current_value(
        Constants.ECCENTRICITIES,
        Constants.ECCENTRICITIES_RATE,
        century_since_epoch
    )

    inclinations = get_current_value(
        Constants.INCLINATIONS,
        Constants.INCLINATIONS_RATE,
        century_since_epoch
    )

    mean_longitudes = get_current_value(
        Constants.MEAN_LONGITUDES,
        Constants.MEAN_LONGITUDES_RATE,
        century_since_epoch
    )

    perihelion_longitudes = get_current_value(
        Constants.PERIHELION_LONGITUDES,
        Constants.PERIHELION_LONGITUDES_RATE,
        century_since_epoch
    )

    ascending_node_longitudes = get_current_value(
        Constants.ASCENDING_NODE_LONGITUDES,
        Constants.ASCENDING_NODE_LONGITUDES_RATE,
        century_since_epoch
    )

    perihelion_arguments = get_perihelion_argument(
        perihelion_longitudes,
        ascending_node_longitudes
    )

    mean_anomalies = get_mean_anomaly(
        mean_longitudes,
        perihelion_longitudes,
        (
            Constants.TERM_B,
            Constants.TERM_C,
            Constants.TERM_S,
            Constants.TERM_F
        ),
        century_since_epoch
    )

    eccentric_anomalies = get_eccentric_anomaly(
        mean_anomalies,
        eccentricities,
        numpy.array(
            [Constants.TOL for num in enumerate(mean_anomalies)],
        )
    )

    heliocentric_coordinates = numpy.column_stack(
        (
            semi_major_axes * (
                cos_wrapper(eccentric_anomalies) - eccentricities
                ),
            semi_major_axes * numpy.sqrt(1 - eccentricities ** 2) * sin_wrapper(eccentric_anomalies),
            numpy.array(
                [0 for num in enumerate(mean_anomalies)],
                numpy.dtype("float64")
            )
        )
    )

    ecliptic_coordinates = get_ecliptic_coordinate(
        heliocentric_coordinates,
        perihelion_arguments,
        ascending_node_longitudes,
        inclinations
    )

    eq_matrix = numpy.matrix(
        [
            [1, 0, 0],
            [
                0,
                cos_wrapper(Constants.OBLIQUITY),
                - sin_wrapper(Constants.OBLIQUITY)
            ],
            [
                0,
                sin_wrapper(Constants.OBLIQUITY),
                cos_wrapper(Constants.OBLIQUITY)
            ]
        ]
    )

    print(numpy.dot(eq_matrix, ecliptic_coordinates.T).T)

if __name__ == "__main__":
    main()
