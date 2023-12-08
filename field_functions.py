from trajectories.trajectory_analysis import *

def SIII_find_axis_unit_vectors(thetas, phis):
    '''
    :parameter: theta - array length N
    :parameter: phi - array length N
    :return: x_hat, y_hat, z_hat - arrays of shape (3,N)
    '''
    x_hats = []
    y_hats = []
    z_hats = []
    for theta, phi in zip(thetas, phis):
        x_hat = np.array([          np.sin(phi),   np.sin(theta) * np.cos(phi),   np.cos(theta) * np.cos(phi)])
        y_hat = np.array([         -np.cos(phi),   np.sin(theta) * np.sin(phi),   np.cos(theta) * np.sin(phi)])
        z_hat = np.array([ np.zeros_like(theta),                -np.cos(theta),                 np.sin(theta)])
        x_hats.append(x_hat)
        y_hats.append(y_hat)
        z_hats.append(z_hat)
    return (x_hats, y_hats, z_hats)

def Galileo_trajectories_SIII_from_CPhiO():
    galileo_callisto_cphio, _ = get_pds_data()
    callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
    galileo_jupiter_SIII = {}
    
    i = 1
    for _orb1, gal_cal_vector in galileo_callisto_cphio.items():
        cal_jup_vector = callisto_jupiter_SIII['orbit%s' % (i)]
        theta_SIII = cal_jup_vector[5]
        phi_SIII = cal_jup_vector[6]
        x_hats_SIII, y_hats_SIII, z_hats_SIII = SIII_find_axis_unit_vectors(theta_SIII, phi_SIII)
        gal_cal_cphio_xyz = np.transpose(gal_cal_vector[1:4])
        x_new = [] ; y_new = [] ; z_new = []
        for gal_cal_cphio_xyz_i, x_hat_SIII, y_hat_SIII, z_hat_SIII in zip(gal_cal_cphio_xyz, x_hats_SIII, y_hats_SIII, z_hats_SIII):
            new_x = np.dot(x_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; x_new.append(new_x)
            new_y = np.dot(y_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; y_new.append(new_y)
            new_z = np.dot(z_hat_SIII, np.transpose(gal_cal_cphio_xyz_i)) ; z_new.append(new_z)
        gal_cal_SIII_cart = np.array([x_new, y_new, z_new])
        gal_jup_SIII_cart = gal_cal_SIII_cart + cal_jup_vector[1:4]
        gal_jup_SIII_spher = cartesian_to_spherical(gal_jup_SIII_cart.transpose())
        gal_jup_SIII_vector = np.c_[cal_jup_vector[0], gal_jup_SIII_cart.transpose()]
        gal_jup_SIII_vector = np.c_[gal_jup_SIII_vector, gal_jup_SIII_spher]
        galileo_jupiter_SIII['orbit%s' % (i)] = gal_jup_SIII_vector.transpose()
        i += 1

    return galileo_jupiter_SIII

galileo_callisto_cphio, _ = get_pds_data()
galileo_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()

callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')

colors = ['r', 'y', 'g']

plt.figure()
for i in range(3):
    plt.plot(galileo_jupiter_SIII['orbit%s' % (i+1)][0] - callisto_jupiter_SIII['orbit%s' % (i+1)][0], color=colors[i])
plt.show()

time_diffs = galileo_jupiter_SIII['orbit%s' % (1)][0] - callisto_jupiter_SIII['orbit%s' % (1)][0]
print(np.linalg.norm(time_diffs))

fig, ax = plt.subplots(2, 3)
for i in range(3):
    ax[0,0].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][1] / R_J, color=colors[i])
    ax[0,0].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][1] / R_J, '--', color=colors[i])
ax[0,0].set_title('x')

for i in range(3):
    ax[0,1].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][2] / R_J, color=colors[i])
    ax[0,1].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][2] / R_J, '--', color=colors[i])
ax[0,1].set_title('y')

for i in range(3):
    ax[0,2].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][3] / R_J, color=colors[i])
    ax[0,2].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][3] / R_J, '--', color=colors[i])
ax[0,2].set_title('z')

for i in range(3):
    ax[1,0].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][4] / R_J, color=colors[i])
    ax[1,0].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][4] / R_J, '--', color=colors[i])
ax[1,0].set_title('r')

for i in range(3):
    ax[1,1].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][5], color=colors[i])
    ax[1,1].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][5], '--', color=colors[i])
ax[1,1].set_title('theta')

for i in range(3):
    ax[1,2].plot(galileo_jupiter_SIII['orbit%s' % (i+1)][6], color=colors[i])
    ax[1,2].plot(callisto_jupiter_SIII['orbit%s' % (i+1)][6], '--', color=colors[i])
ax[1,2].set_title('phi')

plt.show()
        

        

    
