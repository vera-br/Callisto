#

import numpy as np
from trajectories.trajectory_analysis import cartesian_to_spherical


RJ = 71492 * 1e3
orbital_dist = 26.3 #* RJ

# rotation angles for the magnetic dipole from the VIP4 model
theta_VIP4 = np.pi * 9.5 / 180
phi_VIP4 = np.pi * 159.2 / 180

# not sure if this works correctly!!
def spherical_coordinates_transformation(sph_coords, distance):
    
    r = sph_coords[0]
    theta = sph_coords[1]
    phi = sph_coords[2]

    # Radial Distance
    r_prime = np.sqrt(r**2 + distance**2 - 2 * r * distance * np.cos(theta))

    # Polar Angle
    cos_theta_prime = (r * np.cos(theta) - distance) / r_prime
    sin_theta_prime = (r * np.sin(theta)) / (r_prime * np.sin(phi))

    # Avoid floating-point errors for arccos and arcsin
    cos_theta_prime = np.clip(cos_theta_prime, -1, 1)
    sin_theta_prime = np.clip(sin_theta_prime, -1, 1)

    theta_prime = np.arccos(cos_theta_prime)
    theta_prime = np.where(sin_theta_prime < 0, 2 * np.pi - theta_prime, theta_prime)

    # Azimuthal Angle
    tan_phi_prime = (r * np.sin(theta) * np.sin(phi)) / (r * np.sin(theta) * np.cos(phi) - distance)
    phi_prime = np.arctan(tan_phi_prime)

    return np.column_stack((r_prime, theta_prime, phi_prime))


def B_khurana_2(orbit_JSO, orbit_SIII_mag):

    # hinging distance in terms of JSO x coordinate
    x = orbit_JSO[1] / RJ

    # cylindrical radial distance of spacecraft from jupiter
    rho = np.sqrt(orbit_SIII_mag[1]**2 + orbit_SIII_mag[2]**2) / RJ
    
    # azimuthal angle measured from prime meridian at 202°
    phi = orbit_SIII_mag[6]

    #distance of spacecraft from centre of current sheet
    Z = orbit_SIII_mag[3] / RJ

    # radial coordinate
    r = orbit_SIII_mag[4] / RJ

    # angular velocity (synodic)
    omegaJ = 2 *np.pi / 10.18 #hr^-1

    # fit parameters from Khurana1992
    x0 = -33.5 #RJ
    rho0 = 33.2 #RJ
    v0 = 37.4 #RJ hr^-1

    # fit parameters from Common model in Khurana1997
    C1 = 80.3
    C2 = 690.4
    C3 = 101.3
    C4 = -1.7
    a1 = 2.49
    a2 = 1.80
    a3 = 2.64
    r01 = 38.0
    rho02 = 2.14 
    rho03 = 12.5
    D1 = 2.01
    D2 = 13.27
    p = 6.26*1e-3
    q = 0.35

    # distance from current sheet to dipole mag equator (eqn 7 & 8)
    delta = np.pi -  omegaJ * rho0/v0 * np.log(np.cosh(rho/rho0))

    Zcs = rho * np.tan(9.6 * np.pi/180) * ( x0/x * np.tanh(x/x0) * np.cos(phi - delta) - np.cos(phi - np.pi))

    # eqn 25
    dZcs_drho = np.tan(9.6 * np.pi/180) * (np.cosh(x/x0)**-2 * np.cos(phi - delta) - rho * x0/x * np.tanh(x/x0) * np.sin(phi - delta) * omegaJ/v0 * np.tanh(rho/rho0) - np.cos(phi - np.pi))

    dZcs_dphi = - rho * np.tan(9.6 * np.pi/180) * (x0/x * np.tanh(x/x0) * np.sin(phi - delta) - np.sin(phi - np.pi))


    # Euler potentials (eqn 24)
    df_drho = -C1 * np.tanh(r01/r)**a1 * np.log(np.cosh((Z-Zcs)/D1)) + C1 * a1 * r01 * rho**2 * r**-3 * np.tanh(r01/r)**(a1-1) * np.cosh(r01/r)**-2 * np.log(np.cosh((Z-Zcs)/D1)) + C1/D1 * rho * np.tanh(r01/r)**a1 * np.tanh((Z-Zcs)/D1) * dZcs_drho + C2 * rho * np.tanh(rho02/rho)**a2 + C3 * rho * np.tanh(rho03/rho)**a3 + C4 * rho

    df_dphi = C1/D1 * rho * np.tanh(r01/r)**a1 * np.tanh((Z-Zcs)/D1) * dZcs_dphi

    df_dz = C1 * a1 * r01 * rho * Z * r**-3 * np.tanh(r01/r)**(a1-1) * np.cosh(r01/r)**-2 * np.log(np.cosh((Z-Zcs)/D1)) - C1/D1 * rho * np.tanh(r01/r)**a1 * np.tanh((Z-Zcs)/D1)

    dg_drho = p * (1 + q * np.tanh((Z-Zcs)/D2)**2) - 2/D2 * p * q * rho * np.tanh((Z-Zcs)/D2) * (np.cosh((Z-Zcs)/D2))**-2 * dZcs_drho

    dg_dphi = 1 - 2/D2 * p * q * rho * np.tanh((Z-Zcs)/D2) * (np.cosh((Z-Zcs)/D2))**-2 * dZcs_dphi

    dg_dz = 2/D2 * p * q * rho * np.tanh((Z-Zcs)/D2) * (np.cosh((Z-Zcs)/D2))**-2


    # external field in SIII (eqn 19)
    Brho = (df_dphi * dg_dz - df_dz * dg_dphi)/rho

    Bphi = df_dz * dg_drho - df_drho * dg_dz

    Bz = (df_drho * dg_dphi - df_dphi * dg_drho)/rho

    B_cyl_SIII_mag = [Brho, Bphi, Bz]


    # rotation matrix to tilt system about y
    rot_matrix_theta = np.array([[ np.cos(theta_VIP4),   0,   -np.sin(theta_VIP4)],
                                 [                  0,   1,                    0],
                                 [ np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]])
    
    # rotation matrix to set prime meridian 
    rot_matrix_phi =np.array( [[ np.cos(phi_VIP4),   np.sin(phi_VIP4),   0],
                               [-np.sin(phi_VIP4),   np.cos(phi_VIP4),   0],
                               [                0,                  0,   1]])
    

    B_cart_SIII = np.empty((0,3))
    for phi_i, B_cyl_i in zip(phi, np.transpose(B_cyl_SIII_mag)):
        
        # rotate and tilt into RH SIII
        B_cyl_SIII_i = np.dot(rot_matrix_theta, np.dot(rot_matrix_phi, B_cyl_i))

        # transformation matrix to convert from cylindrical to cartesian coords
        rot_matrix_cyl_cart = [[ np.cos(phi_i),  -np.sin(phi_i),   0], 
                               [ np.sin(phi_i),   np.cos(phi_i),   0], 
                               [             0,               0,   1]]

        B_cart_SIII_i = np.dot(rot_matrix_cyl_cart, B_cyl_SIII_i)

        # append cartesian values to array
        B_cart_SIII = np.vstack([B_cart_SIII, B_cart_SIII_i])

    # cartesian to spherical transformation in SIII coords.
    B_spher_SIII = cartesian_to_spherical(B_cart_SIII)    
    B_spher_SIII = np.transpose(B_spher_SIII)
    Br, Btheta, Bphi = B_spher_SIII

    # conversion into CPhiO coord. system
    B_spher_cphio = spherical_coordinates_transformation(B_spher_SIII, orbital_dist)

    Bx = B_spher_cphio[:,2] #phi
    By = B_spher_cphio[:,0] #r
    Bz = -(90 - B_spher_cphio[:,1]* 180/np.pi) # -theta but converted to deg

    return np.array([Bx, By, Bz]).transpose()