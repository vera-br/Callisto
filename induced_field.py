# induced field

import numpy as np
from scipy import constants
import scipy.special as sps
from maths import *

mu0 = constants.mu_0
pi = constants.pi

def R_bessel(r1, k):
    z = r1 * k
    J_min_52 = sps.jv(-5 / 2, z)
    J_32 = sps.jv(3 / 2, z)
    J_12 = sps.jv(1 / 2, z)

    return (z * J_min_52) / (3 * J_32 - z * J_12)


def ae_iphi(sigma, omega, rm, r0, r1):
    k = (1 - 1j) * np.sqrt(mu0 * sigma * omega / 2)
    R = R_bessel(r1, k)

    amp = (r0 / rm) ** 3
    z = r0 * k

    J_52 = sps.jv(5 / 2, z)
    J_min_52 = sps.jv(-5 / 2, z)
    J_12 = sps.jv(1 / 2, z)
    J_min_12 = sps.jv(-1 / 2, z)

    return np.real(amp * (R * J_52 - J_min_52) / (R * J_12 - J_min_12))


def B_induced_finite_conductivity(orbit, B_external, sigma, omega, Rm, R0, R1):
    """
    Calculate the induced magnetic field with finite conductivity
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param sigma: conductivity in S
    :param omega: angular frequency of inducing field in
    :param rm: object radius in m
    :param r0: conducting layer outer radius in m
    :param r1: conducting layer inner radius in m
    :return: time evolution array of Bx, By, Bz in nT
    for more info see (Zimmer et al) https://www.sciencedirect.com/science/article/pii/S001910350096456X
    """
    orbit = orbit.transpose()

    A = ae_iphi(sigma, omega, Rm, R0, R1)

    Bind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        position = vector[1:4] 

        M = -(2 * pi / mu0) * A * B_ext * (Rm**3)

        rmag = np.linalg.norm(position)
        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)


def B_induced_infinite(orbit, B_external, Rm, R0):
    """
    Induced magnetic field for infinite conductivity (superconductor)
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :param B_external: array of external field vectors Bx, By, Bz in nT
    :param Rm: object radius in m
    :param R0: conducting layer outer radius in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    orbit = orbit.transpose()

    A = -((R0 / Rm) ** 3)

    B_ind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        pos = vector[1:4] 

        M = -(2 * pi / mu0) * A * 2 * B_ext * (Rm**3)

        rmag = np.linalg.norm(pos)
        rdotM_r = np.dot(pos, M) * pos

        B_ind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        B_ind_evolution.append(B_ind)

    return np.array(B_ind_evolution)
