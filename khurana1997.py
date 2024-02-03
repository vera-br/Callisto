import numpy as np
from trajectories.trajectory_analysis import cartesian_to_spherical
from field_functions import *

R_J = 71492 * 1e3

def B_sheet_khurana(orbit_JSO, orbit_SIII_mag, orbit_SIII):
    O_JSO = orbit_JSO.copy()
    O_SIII_mag = orbit_SIII_mag.copy()
    O_SIII = orbit_SIII.copy()
    
    x_JSO = O_JSO[1] / R_J

    # coords. in magnetodisc frame
    rho = np.sqrt(O_SIII_mag[1]**2 + O_SIII_mag[2]**2) / R_J
    Z = O_SIII_mag[3] / R_J
    r = O_SIII_mag[4] / R_J

    psi = O_SIII_mag[6]
    theta_mag = O_SIII_mag[5]

    # fit constants from Khurana (1997)
    x0 = -33.5 #* R_J
    rho0 = 33.2 #* R_J
    v0 = 37.4 #* R_J hr^-1
    omega_J = 2 * np.pi / 10.1 # hr

    C1 = 80.3 ; C2 = 690.4 ; C3 = 101.3 ; C4 = -1.7
    a1 = 2.49 ; a2 = 1.80 ; a3 = 2.64
    r01 = 38.0 #* R_J
    rho02 = 2.14 #* R_J
    rho03 = 12.5 #* R_J 
    D1 = 2.01 #* R_J 
    D2 = 13.27 #* R_J
    p = 6.26e-3 ; q = 0.35 
    
    delta = np.pi - omega_J * rho0 / v0 * np.log(np.cosh(rho / rho0))
    
    tan_96 = np.tan(np.pi * 9.6 / 180)
    tanh_x_x0 = np.tanh(x_JSO / x0)

    Zcs = rho * tan_96 * (x0 / x_JSO * tanh_x_x0 * np.cos(psi - delta) - np.cos(psi - np.pi))
    dZcs_drho = tan_96 * ((np.cosh(x_JSO / x0))**-2 * np.cos(psi - delta) \
                - rho * x0 / x_JSO * tanh_x_x0 * np.sin(psi - delta) * omega_J / v0 * np.tanh(rho / rho0) \
                - np.cos(psi - np.pi))
    dZcs_dpsi = -rho * tan_96 * (x0 / x_JSO * tanh_x_x0 * np.sin(psi - delta) - np.sin(psi - np.pi))
    
    ln_coshZZcs_D1 = np.log(np.cosh((Z - Zcs) / D1))
    coshZZcs_D2 = np.cosh((Z - Zcs) / D2)
    tanhZZcs_D1 = np.tanh((Z - Zcs) / D1)
    tanhZZcs_D2 = np.tanh((Z - Zcs) / D2)
    tanh_r01_r = np.tanh(r01 / r)
    sech2_r01_r = (np.cosh(r01 / r))**-2
    
    df_drho = -C1 * tanh_r01_r**a1 * ln_coshZZcs_D1 \
              + C1 * a1 * r01 * rho**2 / r**3 * tanh_r01_r**(a1-1) * sech2_r01_r * ln_coshZZcs_D1 \
              + C1 / D1 * rho * tanh_r01_r**a1 * tanhZZcs_D1 * dZcs_drho \
              + C2 * rho * (np.tanh(rho02 / rho))**a2 \
              + C3 * rho * (np.tanh(rho03 / rho))**a3 \
              + C4 * rho

    df_dpsi = C1 * rho / D1 * tanh_r01_r**a1 * tanhZZcs_D1 * dZcs_dpsi

    df_dz = C1 * a1 * r01 * rho * Z / r**3 * tanh_r01_r**(a1 - 1) * sech2_r01_r * ln_coshZZcs_D1 \
            - C1 * rho / D1 * tanh_r01_r**a1 * tanhZZcs_D1

    dg_drho = p * (1 + q * tanhZZcs_D2**2) \
              - 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2 * dZcs_drho

    dg_dpsi = 1 - 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2 * dZcs_dpsi
    
    dg_dz = 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2
    
    # Cylindrical B Field in magnetodisc frame
    B_rho = (df_dpsi * dg_dz - df_dz * dg_dpsi) / rho

    B_psi = df_dz * dg_drho - df_drho * dg_dz

    B_z = (df_drho * dg_dpsi - df_dpsi * dg_drho) / rho
    
    B_cyl = [B_rho, B_psi, B_z]
    
    # angles from SIII frame for cartesian to spherical transformation
    theta_SIII = O_SIII[5]
    phi_SIII = O_SIII[6]

    B_spher_SIII = []
    for psi_i, B_cyl_i, theta_i, phi_i in zip(psi, np.transpose(B_cyl), theta_SIII, phi_SIII):
        rot_matrix_cyl_cart = [[ np.cos(psi_i),  -np.sin(psi_i),   0], 
                               [ np.sin(psi_i),   np.cos(psi_i),   0], 
                               [             0,               0,   1]]

        # cylindrical to cartesian transformation in magnetodisc coords.
        B_cart_mag_i = np.dot(rot_matrix_cyl_cart, B_cyl_i)

        # rotation of dipole relating to untilting coord. axis by theta_VIP4 then unrotating by phi_VIP4
        B_cart_SIII_i = rotation_SIII_mag_to_SIII(B_cart_mag_i)
        
        rot_matrix_cart_spher = [[ np.cos(phi_i) * np.sin(theta_i),     np.sin(phi_i) * np.sin(theta_i),   np.cos(theta_i)], 
                                 [ np.cos(phi_i) * np.cos(theta_i),     np.sin(phi_i) * np.cos(theta_i),  -np.sin(theta_i)], 
                                 [                  -np.sin(phi_i),                       np.cos(phi_i),                 0]]
        
        # cartesian to spherical transformation in SIII coords.
        B_spher_SIII_i = np.dot(rot_matrix_cart_spher, B_cart_SIII_i)
        B_spher_SIII.append(B_spher_SIII_i)
    
    B_spher_SIII = np.transpose(B_spher_SIII)
    Br, Btheta, Bphi = B_spher_SIII
    # conversion into CPhiO coord. system
    Bx = Bphi
    By = -Br
    Bz = -Btheta
    return np.array([Bx, By, Bz]).transpose()

def B_sheet_khurana2(orbit_JSO, orbit_SIII_mag, model="common"):
    """
    Can choose fit parameters based on "Pioneer10", "Voyager1", "Voyager2", or a combination of all.
    Default is all i.e. "common" model.
    """

    O_JSO = orbit_JSO.copy()
    O_SIII_mag = orbit_SIII_mag.copy()
    
    x_JSO = O_JSO[1] / R_J

    # coords. in magnetodisc frame
    rho = np.sqrt(O_SIII_mag[1]**2 + O_SIII_mag[2]**2) / R_J
    Z = O_SIII_mag[3] / R_J
    r = O_SIII_mag[4] / R_J

    psi = O_SIII_mag[6]
    theta_mag = O_SIII_mag[5]

    # fit constants from Khurana (1997)
    x0 = -33.5 #* R_J
    rho0 = 33.2 #* R_J
    v0 = 37.4 #* R_J hr^-1
    omega_J = 2 * np.pi / 10.1 # hr

    # Fit parameters fom Pioneer 10, Voyager 1 & 2, and all combined
    C1 = [70.2, 100.8, 84.8, 80.3]
    C2 = [1369.9, 916.4, 1034.7, 690.4]
    C3 = [33.4, 70.6, 39.3, 101.3]
    C4 = [-1.1, -1.3, -1.5, -1.7]
    a1 = [3.27, 2.61, 2.04, 2.49]
    a2 = [2.06, 1.62, 1.75, 1.80]
    a3 = [7.55, 9.56, 4.43, 2.64]
    r01 = [44.1, 35.6, 36.7, 38.0] #R_J
    rho02 = [2.55, 1.78, 1.99, 2.14] #R_J
    rho03 = [32.8, 27.9, 19.7, 12.5] #R_J
    D1 = [1.83, 2.16, 2.21, 2.01] #R_J
    D2 = [20.6, 19.23, 16.81, 13.27]
    p = [6.66*1e-3, 7.48*1e-3, 5.73*1e-3, 6.26*1e-3]
    q = [0.32, 0.33, 0.35, 0.35]

    if model == "Pioneer10":
        C1 = C1[0]; C2 = C2[0]; C3 = C3[0]; C4 = C4[0]
        a1 = a1[0]; a2 = a2[0]; a3 = a3[0]
        r01 = r01[0]
        rho02 = rho02[0]; rho03 = rho03[0]
        D1 = D1[0]; D2 = D2[0]
        p = p[0]; q = q[0]

    elif model == "Voyager1":
        C1 = C1[1]; C2 = C2[1]; C3 = C3[1]; C4 = C4[1]
        a1 = a1[1]; a2 = a2[1]; a3 = a3[1]
        r01 = r01[1]
        rho02 = rho02[1]; rho03 = rho03[1]
        D1 = D1[1]; D2 = D2[1]
        p = p[1]; q = q[1]
    
    elif model == "Voyager2":
        C1 = C1[2]; C2 = C2[2]; C3 = C3[2]; C4 = C4[2]
        a1 = a1[2]; a2 = a2[2]; a3 = a3[2]
        r01 = r01[2]
        rho02 = rho02[2]; rho03 = rho03[2]
        D1 = D1[2]; D2 = D2[2]
        p = p[2]; q = q[2]

    else:
        C1 = C1[3]; C2 = C2[3]; C3 = C3[3]; C4 = C4[3]
        a1 = a1[3]; a2 = a2[3]; a3 = a3[3]
        r01 = r01[3]
        rho02 = rho02[3]; rho03 = rho03[3]
        D1 = D1[3]; D2 = D2[3]
        p = p[3]; q = q[3]
    
    delta = np.pi - omega_J * rho0 / v0 * np.log(np.cosh(rho / rho0))
    
    tan_96 = np.tan(np.pi * 9.6 / 180)
    tanh_x_x0 = np.tanh(x_JSO / x0)

    Zcs = rho * tan_96 * (x0 / x_JSO * tanh_x_x0 * np.cos(psi - delta) - np.cos(psi - np.pi))
    dZcs_drho = tan_96 * ((np.cosh(x_JSO / x0))**-2 * np.cos(psi - delta) \
                - rho * x0 / x_JSO * tanh_x_x0 * np.sin(psi - delta) * omega_J / v0 * np.tanh(rho / rho0) \
                - np.cos(psi - np.pi))
    dZcs_dpsi = -rho * tan_96 * (x0 / x_JSO * tanh_x_x0 * np.sin(psi - delta) - np.sin(psi - np.pi))
    
    ln_coshZZcs_D1 = np.log(np.cosh((Z - Zcs) / D1))
    coshZZcs_D2 = np.cosh((Z - Zcs) / D2)
    tanhZZcs_D1 = np.tanh((Z - Zcs) / D1)
    tanhZZcs_D2 = np.tanh((Z - Zcs) / D2)
    tanh_r01_r = np.tanh(r01 / r)
    sech2_r01_r = (np.cosh(r01 / r))**-2
    
    df_drho = -C1 * tanh_r01_r**a1 * ln_coshZZcs_D1 \
              + C1 * a1 * r01 * rho**2 / r**3 * tanh_r01_r**(a1-1) * sech2_r01_r * ln_coshZZcs_D1 \
              + C1 / D1 * rho * tanh_r01_r**a1 * tanhZZcs_D1 * dZcs_drho \
              + C2 * rho * (np.tanh(rho02 / rho))**a2 \
              + C3 * rho * (np.tanh(rho03 / rho))**a3 \
              + C4 * rho

    df_dpsi = C1 * rho / D1 * tanh_r01_r**a1 * tanhZZcs_D1 * dZcs_dpsi

    df_dz = C1 * a1 * r01 * rho * Z / r**3 * tanh_r01_r**(a1 - 1) * sech2_r01_r * ln_coshZZcs_D1 \
            - C1 * rho / D1 * tanh_r01_r**a1 * tanhZZcs_D1

    dg_drho = p * (1 + q * tanhZZcs_D2**2) \
              - 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2 * dZcs_drho

    dg_dpsi = 1 - 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2 * dZcs_dpsi
    
    dg_dz = 2 * p * q / D2 * rho * tanhZZcs_D2 * coshZZcs_D2**-2
    
    # Cylindrical B Field in magnetodisc frame
    B_rho = (df_dpsi * dg_dz - df_dz * dg_dpsi) / rho

    B_psi = df_dz * dg_drho - df_drho * dg_dz

    B_z = (df_drho * dg_dpsi - df_dpsi * dg_drho) / rho
    
    B_cyl = [B_rho, B_psi, B_z]
    
    B_cyl_tilted = []
    for B_cyl_i, theta_i in zip(np.transpose(B_cyl), theta_mag):
        rot_matrix_theta = [[ np.sin(theta_i),   0,   np.cos(theta_i)], 
                            [               0,   1,                 0], 
                            [-np.cos(theta_i),   0,   np.sin(theta_i)]]
        B_cyl_tilted_i = np.dot(rot_matrix_theta, B_cyl_i)
        B_cyl_tilted.append(B_cyl_tilted_i)
    
    Brho, Bpsi, Bz = np.transpose(B_cyl_tilted)
    # conversion into CPhiO coord. system
    Bx = Bpsi
    By = -Brho
    Bz = Bz
    return np.array([Bx, By, Bz]).transpose()