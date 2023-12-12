import numpy as np
import matplotlib.pyplot as plt
from numba import njit

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


def B_current_sheet_rho_variation(Icon, D, rho, z, a, R1):

    Brho = []
    Bz = []

    for i, rho_val in enumerate(rho):
        Brho_sub, Bz_sub = small_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z
                                            , a=R1)

        if rho_val > 5:
            Brho_, Bz_ = large_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z, a=a)
            # Brho.append(Brho_ - Brho_sub)
            # Bz.append(Bz_ - Bz_sub)

            Brho.append(Brho_)
            Bz.append(Bz_)

        elif rho_val <= 5:
            Brho_, Bz_ = small_rho_approx(Icon=Icon, D=D, rho=rho_val, z=z, a=a)
            Brho.append(Brho_)
            Bz.append(Bz_)

    return np.array(Brho), np.array(Bz)

