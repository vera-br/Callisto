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

    return amp * (R * J_52 - J_min_52) / (R * J_12 - J_min_12)


def B_induced_finite_conductivity(pos_vectors, Bext_vectors, sigma, omega, rm, r0, r1):
    """
    Calculate the induced magnetic field with finite conductivity
    :param pos_vectors: positions
    :param Bext_vectors: external field vectors
    :param sigma: conductivity
    :param omega: angular frequency of inducing field
    :param rm: object radius
    :param r0: conducting layer outer radius
    :param r1: conducting layer inner radius
    :return: array of induced magnetic field vectors
    for more info see (Zimmer et al) https://www.sciencedirect.com/science/article/pii/S001910350096456X
    """

    Ae_iphi = ae_iphi(sigma, omega, rm, r0, r1)

    Bind_evolution = []
    for bext, position in zip(Bext_vectors, pos_vectors):

        M = -(2 * pi / mu0) * Ae_iphi * bext * (rm**3)
        rmag = np.linalg.norm(position)

        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)


def B_induced_superconductor(pos_vectors, Bext_vectors, rm, r0):
    """
    Same as B_induced_finite_conductivity() but for a superconductor (i.e., sigma -> infinity)
    :param pos_vectors:
    :param Bext_vectors:
    :param rm: object radius
    :param r0: conducting layer outer radius
    :return:
    """
    A = -((r0 / rm) ** 3)

    Bind_evolution = []
    for bext, position in zip(Bext_vectors, pos_vectors):

        M = -(2 * pi / mu0) * A * bext * (rm**3)
        rmag = np.linalg.norm(position)

        rdotM_r = np.dot(position, M) * position

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)