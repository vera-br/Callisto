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



def B_induced_infinite(orbit, B_external, Rm, R0):
    """
    Induced magnetic field for infinite conductivity (superconductor)
    :param orbit:
    :param B_esxternal:
    :param Rm: object radius
    :param R0: conducting layer outer radius
    :return:
    """
    orbit = orbit.transpose()

    A = -((R0 / Rm) ** 3)

    B_ind_evolution = []
    for B_ext, vector in zip(B_external, orbit):

        pos = vector[1:4] 

        M = -(2 * pi / mu0) * A * B_ext * (Rm**3)

        rmag = np.linalg.norm(pos)
        rdotM_r = np.dot(pos, M) * pos

        B_ind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        B_ind_evolution.append(B_ind)

    return np.array(B_ind_evolution)
