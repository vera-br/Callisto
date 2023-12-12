# plasma sheet

import numpy as np

# matthew and ciaran's files
from maths import cylindrical_to_cartesian

RJ = 71492e3

# parameters for cylindrical plasma sheet
R0 = 7.8  # disc inner radius (RJ)
R1 = 51.4  # disc outer radius (RJ)
D = 3.6  # disc half thickness (RJ)
Icon = 139.6  # current constant = mu0 * I / 2 (nT)
#thetaD = np.radians(9.3)  # disc normal from rotation axis (radians)
#phiD = np.radians(204.2)  # azimuth angle of disc normal (radians)

def B_currents_interior(I_constant, rho, z, D, a_inner, a_outer):
    """
    Calculates the radial and vertical magnetic fields due do currents in the plasma sheet where
    the current profile is I(p) = I0 / p

    The coordinate system is:
        rho: radial position perpendicular to magnetic axis
        phi: azimuthal direction
        z: aligned with magnetic axis

    :param I_constant: current constant (mu0 * I / 2) in nT
    :param rho: radial position in m??
    :param z: vertical position in m??
    :param D: disk half thickness in RJ
    :param a_inner: disk inner edge in RJ
    :param a_outer: disk outer edge in RJ
    :return:
    """

    def inner(I_constant, rho, z, D, a_inner):

        F1 = np.sqrt((z - D) ** 2 + a_inner ** 2)
        F2 = np.sqrt((z + D) ** 2 + a_inner ** 2)

        Brho_term1 = (1 / rho) * (F1 - F2 + (2 * z))
        Brho_term2 = - ((a_inner ** 2 * rho) / 4) * ((1 / (F1 ** 3)) - (1 / (F2 ** 3)))

        Brho = I_constant * (Brho_term1 + Brho_term2)

        Bz_term1 = 2 * D * (1 / np.sqrt(z**2 + rho**2))
        Bz_term2 = (a_inner ** 2 / 4) * (((z - D) / (F1 ** 3)) - ((z + D) / (F2 ** 3)))

        Bz = I_constant * (Bz_term1 - Bz_term2)

        return Brho, Bz

    def outer(I_constant, rho, z, D, a_outer):

        F1 = np.sqrt((z - D) ** 2 + a_outer ** 2)
        F2 = np.sqrt((z + D) ** 2 + a_outer ** 2)

        Brho_term1 = (1 / rho) * (F1 - F2 + (2 * z))
        Brho_term2 = - ((a_outer ** 2 * rho) / 4) * ((1 / (F1 ** 3)) - (1 / (F2 ** 3)))

        Brho = I_constant * (Brho_term1 + Brho_term2)

        Bz_term1 = 2 * D * (1 / np.sqrt(z ** 2 + rho ** 2))
        Bz_term2 = (a_outer ** 2 / 4) * (((z - D) / (F1 ** 3)) - ((z + D) / (F2 ** 3)))

        Bz = I_constant * (Bz_term1 - Bz_term2)

        return Brho, Bz

    Brho_inner, Bz_inner = inner(I_constant, rho, z, D, a_inner)
    Brho_outer, Bz_outer = outer(I_constant, rho, z, D, a_outer)

    Brho = Brho_inner - Brho_outer
    Bz = Bz_inner - Bz_outer

    return Brho, Bz

def B_disk(orbit, R0, R1, D, I_constant):
    """
    :param orbit: juice wrt jupiter in SIII
    :param R0: inner disk radius in RJ
    :param R1: outer disk radius in RJ
    :param D: disk half thickness in RJ
    :param I_constant: current constant (mu0 * I / 2) in nT
    :return:
    """

    orbit = orbit.transpose()

    X = []
    Y = []
    Z = []

    for vector in orbit:
        X.append(vector[1] / RJ)
        Y.append(vector[2] / RJ)
        Z.append(vector[3] / RJ)
        
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    
    rho = np.sqrt(X**2 + Y**2)

    Brho_inner, Bz_inner = B_currents_interior(
        I_constant=I_constant, rho=rho, z=Z, D=D, a_inner=R0, a_outer=R1
    )
    Brho_outer, Bz_outer = B_currents_interior(
        I_constant=I_constant, rho=rho, z=Z, D=D, a_inner=R0, a_outer=R1
    )

    Brho = Brho_inner - Brho_outer
    Bz = Bz_inner - Bz_outer
    Bphi = np.full_like(Brho, 0)
    B = np.array([Brho, Bphi, Bz]).transpose()

    Bcart = []
    for vector in B:
        Bcart.append(cylindrical_to_cartesian(A=vector, theta=-np.pi / 2))
    Bcart = np.array(Bcart)

    return rho, Z, Bcart