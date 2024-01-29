import numpy as np
from trajectories.trajectory_analysis import cartesian_to_spherical

R_J = 71492 * 1e3

# rotation angles for the magnetic dipole from the VIP4 model
theta_VIP4 = -np.pi * 9.5 / 180
phi_VIP4 = -np.pi * 159.2 / 180

def convert_SIII_to_SIII_mag(orbit_SIII):
    _orbit_SIII_mag = orbit_SIII
    rot_matrix_theta = [[ np.cos(theta_VIP4),   0,   np.sin(theta_VIP4)], 
                        [                  0,   1,                    0], 
                        [-np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]]
    
    rot_matrix_phi = [[ np.cos(phi_VIP4),  -np.sin(phi_VIP4),   0], 
                      [ np.sin(phi_VIP4),   np.cos(phi_VIP4),   0], 
                      [                0,                  0,   1]]
    _orbit_SIII_mag[1:4] = np.dot(rot_matrix_theta, np.dot(rot_matrix_phi, _orbit_SIII_mag[1:4]))
    _orbit_SIII_mag[4:] = cartesian_to_spherical(_orbit_SIII_mag[1:4].transpose()).transpose()
    return _orbit_SIII_mag

def convert_SIII_mag_to_SIII(orbit_SIII_mag):
    _orbit_SIII = orbit_SIII_mag
    rot_matrix_theta = [[ np.cos(theta_VIP4),   0,   -np.sin(theta_VIP4)], 
                                     [                  0,   1,                     0], 
                                     [ np.sin(theta_VIP4),   0,    np.cos(theta_VIP4)]]
    rot_matrix_phi = [[ np.cos(phi_VIP4),   np.sin(phi_VIP4),   0], 
                                   [-np.sin(phi_VIP4),   np.cos(phi_VIP4),   0], 
                                   [                0,                  0,   1]]
    _orbit_SIII[1:4] = np.dot(rot_matrix_phi, np.dot(rot_matrix_theta, _orbit_SIII[1:4]))
    _orbit_SIII[4:] = cartesian_to_spherical(_orbit_SIII[1:4].transpose()).transpose()
    return _orbit_SIII


def B_sheet_khurana(orbit_JSO, orbit_SIII_mag, orbit_SIII):
    x_JSO = orbit_JSO[1] / R_J

    # coords. in magnetodisc frame
    rho = np.sqrt(orbit_SIII_mag[1]**2 + orbit_SIII_mag[2]**2) / R_J
    Z = orbit_SIII_mag[3] / R_J
    r = orbit_SIII_mag[4] / R_J

    psi = orbit_SIII_mag[6]

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
    
    rot_matrix_theta = [[ np.cos(theta_VIP4),   0,  -np.sin(theta_VIP4)], 
                        [                  0,   1,                    0], 
                        [ np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]]
    
    rot_matrix_phi = [[ np.cos(phi_VIP4),   np.sin(phi_VIP4),   0], 
                      [-np.sin(phi_VIP4),   np.cos(phi_VIP4),   0], 
                      [                0,                  0,   1]]
    
    # angles from SIII frame for cartesian to spherical transformation
    theta_SIII = orbit_SIII[5]
    phi_SIII = orbit_SIII[6]

    B_spher_SIII = []
    for psi_i, B_cyl_i, theta_i, phi_i in zip(psi, np.transpose(B_cyl), theta_SIII, phi_SIII):
        # rot_matrix_cyl_cart = np.transpose([[ np.cos(psi_i),  -np.sin(psi_i),   0], 
        #                                     [ np.sin(psi_i),   np.cos(psi_i),   0], 
        #                                     [             0,               0,   1]])
        rot_matrix_cyl_cart = [[ np.cos(psi_i),  -np.sin(psi_i),   0], 
                               [ np.sin(psi_i),   np.cos(psi_i),   0], 
                               [             0,               0,   1]]

        # cylindrical to cartesian transformation in magnetodisc coords.
        B_cart_mag_i = np.dot(rot_matrix_cyl_cart, B_cyl_i)

        # rotation of dipole relating to untilting coord. axis by theta_VIP4 then unrotating by phi_VIP4
        B_cart_SIII_i = np.dot(rot_matrix_phi, np.dot(rot_matrix_theta, B_cart_mag_i))
        
        rot_matrix_cart_spher = [[ np.cos(phi_i) * np.sin(theta_i),     np.sin(phi_i) * np.sin(theta_i),   np.cos(theta_i)], 
                                 [ np.cos(phi_i) * np.cos(theta_i),     np.sin(phi_i) * np.cos(theta_i),  -np.sin(theta_i)], 
                                 [                  -np.sin(phi_i),                       np.cos(phi_i),                 0]]
        # rot_matrix_cart_spher = np.transpose([[ np.cos(phi_i) * np.sin(theta_i),     np.sin(phi_i) * np.sin(theta_i),   np.cos(theta_i)], 
        #                                       [ np.cos(phi_i) * np.cos(theta_i),     np.sin(phi_i) * np.cos(theta_i),  -np.sin(theta_i)], 
        #                                       [                  -np.sin(phi_i),                       np.cos(phi_i),                 0]])
        
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


    