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

def ae_iphi_multi_layer(conductivities, r, l, omega):
    def Fs(l, r, k):
        rk = r * k
        sqrt = np.sqrt(pi / (2 * rk))
        I_v = sps.iv(l + 0.5, rk)
        dI_v = sps.ivp(l + 0.5, rk)
        K_v = sps.kv(l + 0.5, rk)
        dK_v = sps.kvp(l + 0.5, rk)
        F1 = sqrt * I_v
        F2 = sqrt * K_v
        dF1 = -0.5 * F1 / r + sqrt * dI_v
        dF2 = -0.5 * F2 / r + sqrt * dK_v
        return F1, F2, dF1, dF2
    
    k = np.sqrt(-1j * omega * mu0 * conductivities)

    k1 = k[0]
    k2 = k[1]
    r1 = r[0]
    r2 = r[1]
    F11, F21, dF11, dF21 = Fs(l, r1, k1)
    F12, F22, dF12, dF22 = Fs(l, r1, k2)
    Djmin1_Cjmin1 = (F12 / F22) * ( ( (dF12 / F12) - (dF11 / F11) ) / ( (dF11 / F11) - (dF22, F22)) )

    for i in range(1, len(r) - 1):
       k_jmin1 = k[i]
       k_j = k[i + 1]
       r_jmin1 = r[i]
       r_j = r[i + 1]
       F11, F21, dF11, dF21 = Fs(l, r_jmin1, k_jmin1)
       F12, F22, dF12, dF22 = Fs(l, r_j, k_j)
       numer = (dF12 / F12) - (dF11 / F11) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF12 / F12) - (dF21 / F21))
       denom = (dF11 / F11) - (dF22 / F22) + Djmin1_Cjmin1 * (F21 / F11) * ( (dF21 / F21) - (dF22 / F22))
       Djmin1_Cjmin1 = (F12 / F22) * (numer / denom)
    
    R = r[-1]
    k_J = k[-1]
    F1J, F2J, dF1J, dF2J = Fs(l, R, k_J)
    numer = (dF1J / F1J) - (l + 1) + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) - (l + 1) )
    denom = (dF1J / F1J) + l + Djmin1_Cjmin1 * (F2J / F1J) * ( (dF2J / F2J) + l )
    
    Ae_iphi = numer / denom

    return Ae_iphi