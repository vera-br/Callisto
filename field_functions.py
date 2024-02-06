from trajectories.trajectory_analysis import *

# def SIII_find_axis_unit_vectors(thetas, phis):
#     '''
#     :parameter: theta - array length N
#     :parameter: phi - array length N
#     :return: x_hat, y_hat, z_hat - arrays of shape (3,N)
#     '''
#     x_hats = []
#     y_hats = []
#     z_hats = []
#     for theta, phi in zip(thetas, phis):
#         x_hat = np.array([          np.sin(phi),   np.sin(theta) * np.cos(phi),   np.cos(theta) * np.cos(phi)])
#         y_hat = np.array([         -np.cos(phi),   np.sin(theta) * np.sin(phi),   np.cos(theta) * np.sin(phi)])
#         z_hat = np.array([ np.zeros_like(theta),                -np.cos(theta),                 np.sin(theta)])
#         x_hats.append(x_hat)
#         y_hats.append(y_hat)
#         z_hats.append(z_hat)
#     return (x_hats, y_hats, z_hats)

# def convert_CPhiO_to_SIII(O_cal_SIII, O_cphio):
#     orbit_cal_SIII = O_cal_SIII.copy()
#     orbit_cphio = O_cphio.copy()

#     theta_SIII = orbit_cal_SIII[5]
#     phi_SIII = orbit_cal_SIII[6]
#     x_hats_SIII, y_hats_SIII, z_hats_SIII = SIII_find_axis_unit_vectors(theta_SIII, phi_SIII)
#     gal_cal_cphio_xyz = np.transpose(orbit_cphio[1:4])
#     x_new = [] ; y_new = [] ; z_new = []
#     for gal_cal_cphio_xyz_i, x_hat_SIII, y_hat_SIII, z_hat_SIII in zip(gal_cal_cphio_xyz, x_hats_SIII, y_hats_SIII, z_hats_SIII):
#         new_x = np.dot(x_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; x_new.append(new_x)
#         new_y = np.dot(y_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; y_new.append(new_y)
#         new_z = np.dot(z_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; z_new.append(new_z)
#     gal_cal_SIII_cart = np.array([x_new, y_new, z_new])
#     gal_jup_SIII_cart = gal_cal_SIII_cart + orbit_cal_SIII[1:4]
#     gal_jup_SIII_spher = cartesian_to_spherical(gal_jup_SIII_cart.transpose())
#     gal_jup_SIII_vector = np.c_[orbit_cal_SIII[0], gal_jup_SIII_cart.transpose()]
#     gal_jup_SIII_vector = np.c_[gal_jup_SIII_vector, gal_jup_SIII_spher]
#     return gal_jup_SIII_vector

def convert_CPhiO_to_SIII_(O_cal_SIII_, O_cphio):
    orbit_cal_SIII_ = O_cal_SIII_.copy()
    orbit_cphio = O_cphio.copy()

    thetas = orbit_cal_SIII_[5]
    phis = orbit_cal_SIII_[6]
    gal_cal_cphio_xyz = np.transpose(orbit_cphio[1:4])
    gal_cal_SIII_cart = []
    for gal_cal_cphio_xyz_i, theta, phi in zip(gal_cal_cphio_xyz, thetas, phis):
        rotation_matrix = [[          np.sin(phi),   np.sin(theta) * np.cos(phi),   np.cos(theta) * np.cos(phi)],
                           [         -np.cos(phi),   np.sin(theta) * np.sin(phi),   np.cos(theta) * np.sin(phi)],
                           [                    0,                -np.cos(theta),                 np.sin(theta)]]
        gal_cal_SIII_cart_i = np.dot(rotation_matrix, gal_cal_cphio_xyz_i)
        gal_cal_SIII_cart.append(gal_cal_SIII_cart_i)
    gal_cal_SIII_cart = np.transpose(gal_cal_SIII_cart)
    gal_jup_SIII_cart = gal_cal_SIII_cart + orbit_cal_SIII_[1:4]
    gal_jup_SIII_spher = cartesian_to_spherical(gal_jup_SIII_cart.transpose())
    gal_jup_SIII_vector = np.c_[orbit_cal_SIII_[0], gal_jup_SIII_cart.transpose()]
    gal_jup_SIII_vector = np.c_[gal_jup_SIII_vector, gal_jup_SIII_spher]
    return gal_jup_SIII_vector

def Galileo_trajectories_SIII_from_CPhiO():
    galileo_callisto_cphio, _ = get_pds_data()
    callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
    galileo_jupiter_SIII = {}
    
    i = 1
    for _orb1, gal_cal_vector in galileo_callisto_cphio.items():
        cal_jup_SIII_vector = callisto_jupiter_SIII['orbit%s' % (i)]
        gal_jup_SIII_vector = convert_CPhiO_to_SIII_(cal_jup_SIII_vector, gal_cal_vector)
        galileo_jupiter_SIII['orbit%s' % (i)] = gal_jup_SIII_vector.transpose()
        i += 1
    return galileo_jupiter_SIII

def Galileo_trajectories_SIII_mag_from_CPhiO():
    galileo_callisto_cphio, _ = get_pds_data()
    callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')
    galileo_jupiter_SIII_mag = {}
    
    i = 1
    for _orb1, gal_cal_vector in galileo_callisto_cphio.items():
        cal_jup_SIII_mag_vector = callisto_jupiter_SIII_mag['orbit%s' % (i)]
        gal_jup_SIII_mag_vector = convert_CPhiO_to_SIII_(cal_jup_SIII_mag_vector, gal_cal_vector)
        galileo_jupiter_SIII_mag['orbit%s' % (i)] = gal_jup_SIII_mag_vector.transpose()
        i += 1
    return galileo_jupiter_SIII_mag

# def Galileo_trajectories_SIII_from_CPhiO():
#     galileo_callisto_cphio, _ = get_pds_data()
#     callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
#     galileo_jupiter_SIII = {}
    
#     i = 1
#     for _orb1, gal_cal_vector in galileo_callisto_cphio.items():
#         cal_jup_vector = callisto_jupiter_SIII['orbit%s' % (i)]
#         theta_SIII = cal_jup_vector[5]
#         phi_SIII = cal_jup_vector[6]
#         x_hats_SIII, y_hats_SIII, z_hats_SIII = SIII_find_axis_unit_vectors(theta_SIII, phi_SIII)
#         gal_cal_cphio_xyz = np.transpose(gal_cal_vector[1:4])
#         x_new = [] ; y_new = [] ; z_new = []
#         for gal_cal_cphio_xyz_i, x_hat_SIII, y_hat_SIII, z_hat_SIII in zip(gal_cal_cphio_xyz, x_hats_SIII, y_hats_SIII, z_hats_SIII):
#             new_x = np.dot(x_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; x_new.append(new_x)
#             new_y = np.dot(y_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; y_new.append(new_y)
#             new_z = np.dot(z_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; z_new.append(new_z)
#         gal_cal_SIII_cart = np.array([x_new, y_new, z_new])
#         gal_jup_SIII_cart = gal_cal_SIII_cart + cal_jup_vector[1:4]
#         gal_jup_SIII_spher = cartesian_to_spherical(gal_jup_SIII_cart.transpose())
#         gal_jup_SIII_vector = np.c_[cal_jup_vector[0], gal_jup_SIII_cart.transpose()]
#         gal_jup_SIII_vector = np.c_[gal_jup_SIII_vector, gal_jup_SIII_spher]
#         galileo_jupiter_SIII['orbit%s' % (i)] = gal_jup_SIII_vector.transpose()
#         i += 1

#     return galileo_jupiter_SIII

# rotation angles for the magnetic dipole from the VIP4 model
theta_VIP4 = -np.pi * 9.5 / 180
phi_VIP4 = -np.pi * 159.2 / 180 #159.2 / 180

def rotation_SIII_to_SIII_mag(vector):
    rot_matrix_theta = [[ np.cos(theta_VIP4),   0,   np.sin(theta_VIP4)], 
                        [                  0,   1,                    0], 
                        [-np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]]
    
    rot_matrix_phi = [[ np.cos(phi_VIP4),  -np.sin(phi_VIP4),   0], 
                      [ np.sin(phi_VIP4),   np.cos(phi_VIP4),   0], 
                      [                0,                  0,   1]]
    rotated_vector = np.dot(rot_matrix_theta, np.dot(rot_matrix_phi, vector))
    return rotated_vector

def rotation_SIII_mag_to_SIII(vector):
    rot_matrix_theta = [[ np.cos(theta_VIP4),   0,  -np.sin(theta_VIP4)], 
                        [                  0,   1,                    0], 
                        [ np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]]
    
    rot_matrix_phi = [[ np.cos(phi_VIP4),   np.sin(phi_VIP4),   0], 
                      [-np.sin(phi_VIP4),   np.cos(phi_VIP4),   0], 
                      [                0,                  0,   1]]
    rotated_vector = np.dot(rot_matrix_phi, np.dot(rot_matrix_theta, vector))
    return rotated_vector

def convert_orbit_SIII_to_SIII_mag(orbit_SIII):
    _orbit_SIII_mag = orbit_SIII.copy()
    _orbit_SIII_mag[1:4] = rotation_SIII_to_SIII_mag(_orbit_SIII_mag[1:4])
    _orbit_SIII_mag[4:] = cartesian_to_spherical(_orbit_SIII_mag[1:4].transpose()).transpose()
    return _orbit_SIII_mag

def convert_orbit_SIII_mag_to_SIII(orbit_SIII_mag):
    _orbit_SIII = orbit_SIII_mag.copy()
    _orbit_SIII[1:4] = rotation_SIII_mag_to_SIII(_orbit_SIII[1:4])
    _orbit_SIII[4:] = cartesian_to_spherical(_orbit_SIII[1:4].transpose()).transpose()
    return _orbit_SIII

def convert_B_to_PDS_CPhiO(B_CPhiO, orbit_cal_jup_SIII, orbit_cal_jup_SIII_CA):
    theta_CA = orbit_cal_jup_SIII_CA[5]
    phi_CA = orbit_cal_jup_SIII_CA[6]
    
    thetas = orbit_cal_jup_SIII[5]
    phis = orbit_cal_jup_SIII[6]
    
    delta_theta = thetas - np.ones_like(thetas) * theta_CA
    delta_phi = phis - np.ones_like(phis) * phi_CA

    print(max(delta_theta))
    print(max(delta_phi))

    B_rot = []
    for B, theta_i, phi_i in zip(B_CPhiO, delta_theta, delta_phi):
        rot_matrix_theta = [[  1,                0,                0],
                            [  0,  np.cos(theta_i), -np.sin(theta_i)],
                            [  0,  np.sin(theta_i),  np.cos(theta_i)]]
        rot_matrix_phi = [[  np.cos(phi_i), -np.sin(phi_i),  0],
                          [  np.sin(phi_i),  np.cos(phi_i),  0],
                          [              0,              0,  1]]
        B_i = np.dot(rot_matrix_phi, np.dot(rot_matrix_theta, B))
        B_rot.append(B_i)
    return np.array(B_rot)
        

        

    
