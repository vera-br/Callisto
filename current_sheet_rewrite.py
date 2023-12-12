import numpy as np
from scipy import constants
from maths import cylindrical_to_cartesian
#from maths import angle_between
#from scipy.spatial.transform import Rotation as ROT
#from current_sheet_edwards import B_current_sheet_static
from maths import general_rotation

# constants
pi = constants.pi
RJ = 71492e3  # Jupiter radius
Cal_Jup_sep = 1.880E9

def B_currents_interior(I_constant, rho, z, D, a_inner, a_outer):
    """
    Calculates the radial and vertical magnetic fields due do currents in the plasma sheet where
    the current profile is I(p) = I0 / p

    The coordinate system is:
        rho: radial position perpendicular to magnetic axis
        phi: azimuthal direction
        z: aligned with magnetic axis

    :param I_constant: current constant (mu0 * I / 2)
    :param rho: radial position
    :param z: vertical position
    :param D: disk half thickness
    :param a_inner: disk inner edge
    :param a_outer: disk outer edge
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

def GphiO_to_JupMag(positions_wrt_ganymede, thetaD):
    """
    Transforms from GphiO to JupMag
    :param positions_wrt_ganymede: shape should be (N, 3)
    :param thetaD: radians
    :return:
    """

    # shift to jupiter centered system
    pos_wrt_Gan = positions_wrt_ganymede + np.array([0, Cal_Jup_sep, 0])

    # swap x and y, y -> rho
    pos_wrt_Gan[:, [0, 1]] = pos_wrt_Gan[:, [1, 0]]
    positions_wrt_jupiter = pos_wrt_Gan

    # rotate about [0, 1, 0]
    positions_wrt_jupiter_mag = general_rotation(
        vector=positions_wrt_jupiter,
        rotation_angle=thetaD,
        axis_of_rotation=np.array([0, 1, 0])
    )

    # covert to units of RJ
    positions_wrt_jupiter_mag /= RJ

    return positions_wrt_jupiter_mag


def small_rho_approx(Icon, D, z, rho, a):
    log_term = (z + D + np.sqrt((z + D) ** 2 + a ** 2)) / (
            z - D + np.sqrt((z - D) ** 2 + a ** 2))
    Bz_term1 = np.log(log_term)

    Bz_term2_1 = (z + D) / (((z + D) ** 2 + a ** 2) ** 1.5)
    Bz_term2_2 = - (z - D) / (((z - D) ** 2 + a ** 2) ** 1.5)

    Bz_term2 = (rho ** 2 / 4) * (Bz_term2_1 + Bz_term2_2)

    Bz = Icon * (Bz_term1 + Bz_term2)

    Brho_term1 = (rho / 2) * ((1 /
                               np.sqrt((z - D) ** 2 + a ** 2)) - (1 / np.sqrt((z + D) ** 2 + a ** 2)))

    Brho_term2_1 = (a ** 2 - 2 * (z - D) ** 2) / ((a ** 2 + (z - D) ** 2) ** 2.5)
    Brho_term2_2 = - (a ** 2 - 2 * (z + D) ** 2) / ((a ** 2 + (z + D) ** 2) ** 2.5)

    Brho_term2 = (rho ** 3 / 16) * (Brho_term2_1 + Brho_term2_2)

    Brho = Icon * (Brho_term1 + Brho_term2)

    return Brho, Bz


def large_rho_approx(Icon, D, z, rho, a):
    log_term = (z + D + np.sqrt((z + D) ** 2 + rho ** 2)) / (
            z - D + np.sqrt((z - D) ** 2 + rho ** 2))
    Bz_term1 = np.log(log_term)

    Bz_term2_1 = (z + D) / (((z + D) ** 2 + rho ** 2) ** 1.5)
    Bz_term2_2 = - (z - D) / (((z - D) ** 2 + rho ** 2) ** 1.5)

    Bz_term2 = (a ** 2 / 4) * (Bz_term2_1 + Bz_term2_2)

    Bz = (Icon * (Bz_term1 + Bz_term2))

    if np.abs(z) <= D:
        Brho_term1 = (1 / rho) * (
                np.sqrt((z - D) ** 2 + rho ** 2) - np.sqrt((z + D) ** 2 + rho ** 2))

        Brho_term2_1 = 1 / (((z + D) ** 2 + rho ** 2) ** 1.5)
        Brho_term2_2 = - 1 / (((z - D) ** 2 + rho ** 2) ** 1.5)

        Brho_term2 = (rho * a ** 2 / 4) * (Brho_term2_1 + Brho_term2_2)

        Brho_term3 = (2 * z) / rho

        Brho = Icon * (Brho_term1 + Brho_term2 + Brho_term3)

    # above the sheet
    elif z >= D:
        Brho_term1_above = (1 / rho) * (
                np.sqrt((z - D) ** 2 + rho ** 2) - np.sqrt((z + D) ** 2 + rho ** 2))

        Brho_term2_1_above = 1 / (((z + D) ** 2 + rho ** 2) ** 1.5)
        Brho_term2_2_above = - 1 / (((z - D) ** 2 + rho ** 2) ** 1.5)

        Brho_term2_above = (rho * a ** 2 / 4) * (Brho_term2_1_above + Brho_term2_2_above)

        Brho_term3_above = (2 * D) / rho

        Brho = Icon * (Brho_term1_above + Brho_term2_above + Brho_term3_above)

    # below the sheet
    elif z < - D:
        Brho_term1_below = (1 / rho) * (
                np.sqrt((z - D) ** 2 + rho ** 2) - np.sqrt((z + D) ** 2 + rho ** 2))

        Brho_term2_1_below = 1 / (((z + D) ** 2 + rho ** 2) ** 1.5)
        Brho_term2_2_below = - 1 / (((z - D) ** 2 + rho ** 2) ** 1.5)

        Brho_term2_below = (rho * a ** 2 / 4) * (Brho_term2_1_below + Brho_term2_2_below)

        Brho_term3_below = - (2 * D) / rho

        Brho = Icon * (Brho_term1_below + Brho_term2_below + Brho_term3_below)

    return Brho, Bz


def B_current_sheet_static(Icon, D, positions, a, R1):
    Bsheet = []
    for pos in positions:
        rho, _, z = pos

        Brho_sub, Bz_sub = small_rho_approx(Icon=Icon, D=D, rho=rho, z=z, a=R1)
        B_subtract = np.array([Brho_sub, 0, Bz_sub])

        if rho > 5:
            Brho, Bz = large_rho_approx(Icon=Icon, D=D, rho=rho, z=z, a=a)
            B = np.array([Brho, 0, Bz])
            Bsheet.append(B - B_subtract)

        elif rho <= 5:
            Brho, Bz = small_rho_approx(Icon=Icon, D=D, rho=rho, z=z, a=a)
            B = np.array([Brho, 0, Bz])
            Bsheet.append(B - B_subtract)

    return np.array(Bsheet)


def B_current_sheet_mesh(Icon, D, rho, z, a, R1):
    N = len(rho)

    Brho = np.zeros((N, N))
    Bz = np.zeros((N, N))

    for i, rho_val in enumerate(rho):
        for j, z_val in enumerate(z):

            Brho_sub, Bz_sub = small_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z_val
                                                , a=R1)

            if rho_val > 5:
                Brho_, Bz_ = large_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z_val, a=a)
                Brho[i, j] = Brho_ - Brho_sub
                Bz[i, j] = Bz_ - Bz_sub

            elif rho_val <= 5:
                Brho_, Bz_ = small_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z_val, a=a)
                Brho[i, j] = Brho_ - Brho_sub

                Bz[i, j] = Bz_ - Bz_sub

    return Brho, Bz

def B_sheet_mag(positions, R0, R1, D, I_constant):

    X, Y, Z = positions
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
        Bcart.append(cylindrical_to_cartesian(A=vector, theta=-pi / 2))
    Bcart = np.array(Bcart)

    return rho, Z, Bcart