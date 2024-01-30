# plasma sheet

import numpy as np
from scipy import constants
import JupiterMag as jm

# matthew and ciaran's files
from maths import cylindrical_to_cartesian

RJ = 71492e3
pi = constants.pi

# parameters for cylindrical plasma sheet
R0 = 7.8  # disc inner radius (RJ)
R1 = 51.4  # disc outer radius (RJ)
D = 3.6  # disc half thickness (RJ)
Icon = 139.6  # current constant = mu0 * I / 2 (nT)
#thetaD = np.radians(9.3)  # disc normal from rotation axis (radians)
#phiD = np.radians(204.2)  # azimuth angle of disc normal (radians)

def B_currents_interior(I_constant, rho, z, D, a_inner, a_outer, azimuthal_field=False, I_rho=12):
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
        Brho_term3 = []
        for rho_i, z_i in zip(rho, z):
            # in sheet
            if np.abs(z_i) < D:
                Brho_term3_i = (2 * z_i) / rho_i

            # above the sheet
            elif z_i >= D:
                Brho_term3_i = (2 * D) / rho_i

            # below the sheet
            elif z_i <= - D:
                Brho_term3_i = - (2 * D) / rho_i
            
            Brho_term3.append(Brho_term3_i)

        Brho = I_constant * (Brho_term1 + Brho_term2 + Brho_term3)

        Bz_term1 = 2 * D * (1 / np.sqrt(z**2 + rho**2))
        Bz_term2 = (a_inner ** 2 / 4) * (((z - D) / (F1 ** 3)) - ((z + D) / (F2 ** 3)))

        Bz = I_constant * (Bz_term1 - Bz_term2)

        return Brho, Bz

    def outer(I_constant, rho, z, D, a_outer):

        F1 = np.sqrt((z - D) ** 2 + a_outer ** 2)
        F2 = np.sqrt((z + D) ** 2 + a_outer ** 2)

        Brho_term1 = (1 / rho) * (F1 - F2 + (2 * z))
        Brho_term2 = - ((a_outer ** 2 * rho) / 4) * ((1 / (F1 ** 3)) - (1 / (F2 ** 3)))

        # Brho term 3 - dependent on position relative to sheet edges
        Brho_term3 = []
        for rho_i, z_i in zip(rho, z):
            # in sheet
            if np.abs(z_i) < D:
                Brho_term3_i = (2 * z_i) / rho_i

            # above the sheet
            elif z_i >= D:
                Brho_term3_i = (2 * D) / rho_i

            # below the sheet
            elif z_i <= - D:
                Brho_term3_i = - (2 * D) / rho_i
            
            Brho_term3.append(Brho_term3_i)

        Brho = I_constant * (Brho_term1 + Brho_term2 + Brho_term3)

        Bz_term1 = 2 * D * (1 / np.sqrt(z ** 2 + rho ** 2))
        Bz_term2 = (a_outer ** 2 / 4) * (((z - D) / (F1 ** 3)) - ((z + D) / (F2 ** 3)))



        Bz = I_constant * (Bz_term1 - Bz_term2)

        return Brho, Bz
    
    def B_current_sheet_azimuthal(I_rho, D, z, rho):
        '''
        :param I_rho: - mega-amps [MA]
        :param rho: R_J
        :return B_phi: - nano-Teslas [nT] 
        '''
        Bphi = []
        for z_i, rho_i in zip(z, rho): 
            Bphi_term1 = -2.79752 * (I_rho / rho_i)
            if rho_i == 0:
                Bphi_term2 = 0
            elif np.abs(z_i) >= D and rho_i > 0:
                Bphi_term2 = z_i / np.abs(z_i)
            elif np.abs(z_i) < D and rho_i > 0:
                Bphi_term2 = z_i / D
            Bphi_i = Bphi_term1 * Bphi_term2
            Bphi.append(Bphi_i)
        return Bphi

    Brho_inner, Bz_inner = inner(I_constant, rho, z, D, a_inner)
    Brho_outer, Bz_outer = outer(I_constant, rho, z, D, a_outer)
    
    Brho = Brho_inner - Brho_outer
    Bz = Bz_inner - Bz_outer

    if azimuthal_field == False:
        return Brho, Bz
    
    elif azimuthal_field == True:
        Bphi = B_current_sheet_azimuthal(I_rho, D, z, rho)

        return Brho, Bz, Bphi

    

    

def B_disk(orbit, R0, R1, D, I_constant, azimuthal_field=False, I_rho=12):
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
    if azimuthal_field == False:
        Brho, Bz = B_currents_interior(I_constant=I_constant, rho=rho, z=Z, D=D, a_inner=R0, a_outer=R1, azimuthal_field=azimuthal_field, I_rho=I_rho)
        Bphi = np.full_like(Brho, 0)
        B = np.array([Brho, Bphi, Bz]).transpose()
    elif azimuthal_field == True:
        Brho, Bz, Bphi = B_currents_interior(I_constant=I_constant, rho=rho, z=Z, D=D, a_inner=R0, a_outer=R1, azimuthal_field=azimuthal_field, I_rho=I_rho)
        B = np.array([Brho, Bphi, Bz]).transpose()

    Bcart = []
    for vector in B:
        Bcart.append(cylindrical_to_cartesian(A=vector, theta=-pi / 2))
    Bcart = np.array(Bcart)

    return rho, Z, Bcart

def B_sheet_Community(orbit_SIII):
    O_SIII = orbit_SIII.copy()
    
    # jm.Con2020.Config(equation_type='analytic', CartesianIn=True, CartesianOut=False)
    # x = orbit_SIII[1] / RJ
    # y = orbit_SIII[2] / RJ
    # z = orbit_SIII[3] / RJ
    # Br, Btheta, Bphi = jm.Con2020.Field(x, y, z)
    
    jm.Con2020.Config(equation_type='analytic', CartesianIn=False, CartesianOut=False)
    r = O_SIII[4] / RJ
    Br, Btheta, Bphi = jm.Con2020.Field(r, O_SIII[5], O_SIII[6])
    
    Bx = Bphi
    By = -Br
    Bz = -Btheta
    return np.array([Bx, By, Bz]).transpose()