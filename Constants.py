"""
Constants.py: Provides constants used by Solaris6K
"""

import numpy

# Semi-major axes and the coresponding rates
# [au, au/century]

SEMI_MAJOR_AXES = numpy.array(
    [
        0.38709843,
        0.72332102,
        1.00000018,
        1.52371234,
        5.20248019,
        9.54149883,
        19.18797948,
        30.06952752,
        39.48686035
    ], numpy.dtype("float64")
)

SEMI_MAJOR_AXES_RATE = numpy.array(
    [
        0.00000000,
        -0.00000026,
        -0.00000003,
        0.00000097,
        -0.00002864,
        -0.00003065,
        -0.00020455,
        0.00006447,
        0.00449751
    ], numpy.dtype("float64")
)

# Eccentricities and the coresponding rates
# [, /century]

ECCENTRICITIES = numpy.array(
    [
        0.20563661,
        0.00676399,
        0.01673163,
        0.09336511,
        0.04853590,
        0.05550825,
        0.04685740,
        0.00895439,
        0.24885238
    ], numpy.dtype("float64")
)

ECCENTRICITIES_RATE = numpy.array(
    [
        0.00002123,
        -0.00005107,
        -0.00003661,
        0.00009149,
        -0.00018026,
        -0.00032044,
        -0.00001550,
        0.00000818,
        0.00006026
    ], numpy.dtype("float64")
)

# Inclinations and the coresponding rates
# [degrees, degrees/century]

INCLINATIONS = numpy.array(
    [
        7.00559432,
        3.39777545,
        -0.00054346,
        1.85181869,
        1.29861416,
        2.49424102,
        0.77298127,
        1.77005520,
        17.14104260
    ], numpy.dtype("float64")
)

INCLINATIONS_RATE = numpy.array(
    [
        -0.00590158,
        0.00043494,
        -0.01337178,
        -0.00724757,
        -0.00322699,
        0.00451969,
        -0.00180155,
        0.00022400,
        0.00000501
    ], numpy.dtype("float64")
)

# Mean longitude and the coresponding rates
# [degrees, degrees/century]

MEAN_LONGITUDES = numpy.array(
    [
        252.25166724,
        181.97970850,
        100.46691572,
        -4.56813164,
        34.33479152,
        50.07571329,
        314.20276625,
        304.22289287,
        238.96535011
    ], numpy.dtype("float64")
)

MEAN_LONGITUDES_RATE = numpy.array(
    [
        149472.67486623,
        58517.81560260,
        35999.37306329,
        19140.29934243,
        3034.90371757,
        1222.11494724,
        428.49512595,
        218.46515314,
        145.18042903
    ], numpy.dtype("float64")
)

# Perihelion longitudes and the coresponding rates
# [degrees, degrees/century]

PERIHELION_LONGITUDES = numpy.array(
    [
        77.45771895,
        131.76755713,
        102.93005885,
        -23.91744784,
        14.27495244,
        92.86136063,
        172.43404441,
        46.68158724,
        224.09702598
    ], numpy.dtype("float64")
)

PERIHELION_LONGITUDES_RATE = numpy.array(
    [
        0.15940013,
        0.05679648,
        0.31795260,
        0.45223625,
        0.18199196,
        0.54179478,
        0.09266985,
        0.01009938,
        -0.00968827
    ], numpy.dtype("float64")
)

# Ascending node longitudes and the coresponding rates
# [degrees, degrees/century]

ASCENDING_NODE_LONGITUDES = numpy.array(
    [
        48.33961819,
        76.67261496,
        -5.11260389,
        49.71320984,
        100.29282654,
        113.63998702,
        73.96250215,
        131.78635853,
        110.30167986
    ], numpy.dtype("float64")
)

ASCENDING_NODE_LONGITUDES_RATE = numpy.array(
    [
        -0.12214182,
        -0.27274174,
        -0.24123856,
        -0.26852431,
        0.13024619,
        -0.25015002,
        0.05739699,
        -0.00606302,
        -0.00809981
    ], numpy.dtype("float64")
)

# Additional terms for mean anomaly calculation for Jupiter through Pluto

TERM_B = numpy.array(
    [
        0,
        0,
        0,
        0,
        -0.00012452,
        0.00025899,
        0.00058331,
        -0.00041348,
        -0.01262724
    ], numpy.dtype("float64")
)

TERM_C = numpy.array(
    [
        0,
        0,
        0,
        0,
        0.06064060,
        -0.13434469,
        -0.97731848,
        0.68346318,
        0
    ], numpy.dtype("float64")
)

TERM_S = numpy.array(
    [
        0,
        0,
        0,
        0,
        -0.35635438,
        0.87320147,
        0.17689245,
        -0.10162547,
        0
    ], numpy.dtype("float64")
)

TERM_F = numpy.array(
    [
        0,
        0,
        0,
        0,
        38.35125000,
        38.35125000,
        7.67025000,
        7.67025000,
        0
    ], numpy.dtype("float64")
)

# Other constants

# J2000.0 in Juilian date
J2000_0 = numpy.float64(2451545.0)
# Day to century
CONVERSION_FACTOR = numpy.float64(36525)
#
ECCENTRICITIY_FACTOR = numpy.float64(180) / numpy.pi
# Kepler's Equation error threshold
TOL = numpy.float64("10E-8")
# J2000 frame obliquity
OBLIQUITY = numpy.float64(23.43928)
