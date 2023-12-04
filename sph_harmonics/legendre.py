import numpy as np

# calculates the Schmidt quasi-normalised Legendre functions for a given theta and returns results in a dictionary

"""
The Schmidt quasi-normalised Legendre functions
See Connerney (1993) - Magnetic Fields of the Outer Planets
:return:
"""


def P10(theta):
    return np.cos(theta)


def P11(theta):
    return np.sin(theta)


def P20(theta):
    return 1.5 * ((np.cos(theta) ** 2) - (1 / 3))


def P21(theta):
    return (3 ** 0.5) * np.cos(theta) * np.sin(theta)


def P22(theta):
    return ((3 ** 0.5) / 2) * (np.sin(theta) ** 2)


def P30(theta):
    return 2.5 * np.cos(theta) * ((np.cos(theta) ** 2) - (9 / 15))


def P31(theta):
    return ((5 * (3 ** 0.5)) / (2 * (2 ** 0.5))) * np.sin(theta) * ((np.cos(theta) ** 2) - (3 / 15))


def P32(theta):
    return (((3 ** 0.5) * (5 ** 0.5)) / 2) * np.cos(theta) * (np.sin(theta) ** 2)


def P33(theta):
    return ((5 ** 0.5) / (2 * 2 ** 0.5)) * np.sin(theta) ** 3


def P40(theta):
    return (35 / 8) * ((np.cos(theta) ** 4) - (30 / 35) * (np.cos(theta) ** 2) + (3 / 35))


def P41(theta):
    return ((7 * 5 ** 0.5) / (2 * 2 ** 0.5)) * np.cos(theta) * np.sin(theta) * ((np.cos(theta) ** 2) - (3 / 7))


def P42(theta):
    return ((7 * 5 ** 0.5) / 4) * (np.sin(theta) ** 2) * ((np.cos(theta) ** 2) - (1 / 7))


def P43(theta):
    return ((7 ** 0.5 * 5 ** 0.5) / (2 * 2 ** 0.5)) * np.cos(theta) * (np.sin(theta) ** 3)


def P44(theta):
    return ((7 ** 0.5 * 5 ** 0.5) / 8) * (np.sin(theta) ** 4)


"""
The derivatives of the Schmidt quasi-normalised Legendre functions
:return:
"""


def dP10_dtheta(theta):
    return - np.sin(theta)


def dP11_dtheta(theta):
    return np.cos(theta)


def dP20_dtheta(theta):
    return - 3 * np.cos(theta) * np.sin(theta)


def dP21_dtheta(theta):
    return - 3 ** 0.5 * ((np.sin(theta) ** 2) - (np.cos(theta) ** 2))


def dP22_dtheta(theta):
    return 3 ** 0.5 * np.cos(theta) * np.sin(theta)


def dP30_dtheta(theta):
    return - 0.5 * np.sin(theta) * (15 * (np.cos(theta) ** 2) - 3)


def dP31_dtheta(theta):
    return - ((3 ** 0.5) / (2 ** 1.5)) * np.cos(theta) * (10 * (np.sin(theta) ** 2) - 5 * (np.cos(theta) ** 2) + 1)


def dP32_dtheta(theta):
    return - ((3 ** 0.5 * 5 ** 0.5) / 2) * np.sin(theta) * ((np.sin(theta) ** 2) - 2 * (np.cos(theta) ** 2))


def dP33_dtheta(theta):
    return ((3 * 5 ** 0.5) / (2 ** 1.5)) * np.cos(theta) * (np.sin(theta) ** 2)


def dP40_dtheta(theta):
    return - 0.5 * np.sin(theta) * (35 * (np.cos(theta) ** 3) - 15 * np.cos(theta))


def dP41_dtheta(theta):
    return - ((5 ** 0.5) / (2 ** 1.5)) * ((np.sin(theta) ** 2) * (21 * (np.cos(theta) ** 2) - 3) -
                                          7 * (np.cos(theta) ** 4) + 3 * (np.cos(theta) ** 2))


def dP42_dtheta(theta):
    return - ((5 ** 0.5 * np.cos(theta) * np.sin(theta)) / 2) * (7 * (np.sin(theta) ** 2) -
                                                                 7 * (np.cos(theta) ** 2) + 1)


def dP43_dtheta(theta):
    return - ((5 ** 0.5 * 7 ** 0.5) / (2 ** 1.5)) * (np.sin(theta) ** 2) * ((np.sin(theta) ** 2) -
                                                                            3 * (np.cos(theta) ** 2))


def dP44_dtheta(theta):
    return ((5 ** 0.5 * 7 ** 0.5) / 2) * np.cos(theta) * (np.sin(theta) ** 3)


def dipole(theta):
    p10 = P10(theta)
    p11 = P11(theta)

    dP10 = dP10_dtheta(theta)
    dP11 = dP11_dtheta(theta)

    leg = {
        "P10": p10,
        "P11": p11,
    }

    # Legendre derivatives
    dleg = {
        "dP10": dP10,
        "dP11": dP11,
    }

    return leg, dleg


def quadrupole(theta):
    p10 = P10(theta)
    p11 = P11(theta)
    p20 = P20(theta)
    p21 = P21(theta)
    p22 = P22(theta)

    dP10 = dP10_dtheta(theta)
    dP11 = dP11_dtheta(theta)
    dP20 = dP20_dtheta(theta)
    dP21 = dP21_dtheta(theta)
    dP22 = dP22_dtheta(theta)

    # Legendre polynomials
    leg = {
        "P10": p10,
        "P11": p11,
        "P20": p20,
        "P21": p21,
        "P22": p22
    }

    # Legendre derivatives
    dleg = {
        "dP10": dP10,
        "dP11": dP11,
        "dP20": dP20,
        "dP21": dP21,
        "dP22": dP22
    }

    return leg, dleg


def octopole(theta):
    p10 = P10(theta)
    p11 = P11(theta)
    p20 = P20(theta)
    p21 = P21(theta)
    p22 = P22(theta)
    p30 = P30(theta)
    p31 = P31(theta)
    p32 = P32(theta)
    p33 = P33(theta)
    p40 = P40(theta)
    p41 = P41(theta)
    p42 = P42(theta)
    p43 = P43(theta)
    p44 = P44(theta)

    dP10 = dP10_dtheta(theta)
    dP11 = dP11_dtheta(theta)
    dP20 = dP20_dtheta(theta)
    dP21 = dP21_dtheta(theta)
    dP22 = dP22_dtheta(theta)
    dP30 = dP30_dtheta(theta)
    dP31 = dP31_dtheta(theta)
    dP32 = dP32_dtheta(theta)
    dP33 = dP33_dtheta(theta)
    dP40 = dP40_dtheta(theta)
    dP41 = dP41_dtheta(theta)
    dP42 = dP42_dtheta(theta)
    dP43 = dP43_dtheta(theta)
    dP44 = dP44_dtheta(theta)

    leg = {
        "P10": p10,
        "P11": p11,
        "P20": p20,
        "P21": p21,
        "P22": p22,
        "P30": p30,
        "P31": p31,
        "P32": p32,
        "P33": p33,
        "P40": p40,
        "P41": p41,
        "P42": p42,
        "P43": p43,
        "P44": p44
    }

    dleg = {
        "dP10": dP10,
        "dP11": dP11,
        "dP20": dP20,
        "dP21": dP21,
        "dP22": dP22,
        "dP30": dP30,
        "dP31": dP31,
        "dP32": dP32,
        "dP33": dP33,
        "dP40": dP40,
        "dP41": dP41,
        "dP42": dP42,
        "dP43": dP43,
        "dP44": dP44
    }

    return leg, dleg
