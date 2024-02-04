from trajectories.trajectory_analysis import *
from field_functions import *

#---------constants-----------


prime_meridian = np.radians(249.2)


#---------functions-----------


def transform_SIII_mag_to_SIII(vector_SIII_mag):
    """
    Transform a vector from SIII Mag to SIII.

    Parameters:
    - vector_SIII_mag: 1D NumPy array representing the vector in SIII Mag coordinates [X_mag, Y_mag, Z_mag].

    Returns:
    - vector_SIII: 1D NumPy array representing the vector in SIII coordinates [X, Y, Z].
    """

    # Transformation matrix from SIII Mag to SIII
    transformation_matrix = np.array([
        [np.sin(prime_meridian), 0, -np.cos(prime_meridian)],
        [0, 1, 0],
        [np.cos(prime_meridian), 0, np.sin(prime_meridian)]
    ])

    # Apply the transformation
    vector_SIII = np.dot(transformation_matrix, vector_SIII_mag)

    return vector_SIII



#---------load data-----------


juice_wrt_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_jupiter_SIII = get_spice_data('juice', 'jupiter', 'SIII', 'J')
juice_jupiter_SIII_mag = get_spice_data('juice', 'jupiter', 'SIII_mag', 'J')

flyby_n = 2

orbit_cphio = juice_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = juice_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = juice_jupiter_SIII_mag["orbit%s" % (flyby_n)]

phi_SIII_mag = orbit_SIII_mag[6]


#---------coordinate transformation-----------

# rotation matrix to tilt system about y
rot_matrix_theta = np.array([[ np.cos(theta_VIP4),   0,   np.sin(theta_VIP4)],
                                [                  0,   1,                    0],
                                [-np.sin(theta_VIP4),   0,   np.cos(theta_VIP4)]])

# rotation matrix to set prime meridian 
rot_matrix_phi =np.array( [[ np.cos(phi_VIP4),  -np.sin(phi_VIP4),   0],
                            [ np.sin(phi_VIP4),   np.cos(phi_VIP4),   0],
                            [                0,                  0,   1]])

calculated_SIII = np.empty((0,3))
for phi_i, SIII_mag_i in zip(phi_SIII_mag, np.transpose(orbit_SIII_mag[1:4])):
    # rotate and tilt into RH SIII
    orbit_SIII_i = np.dot(rot_matrix_phi, np.dot(rot_matrix_theta, SIII_mag_i))

    calculated_SIII = np.vstack([calculated_SIII, orbit_SIII_i])

# orbit_SIII = [rotation_SIII_mag_to_SIII(vector) for vector in orbit_SIII_mag[:, 1:4].T]

# x = orbit_SIII[6] #phi
# y = -orbit_SIII[4] #r
# z = -orbit_SIII[5] #theta

# t = orbit_SIII_mag[0].transpose()
# print(np.shape(t))
# print(np.shape(x))

transformed = np.column_stack([orbit_SIII_mag[0], calculated_SIII]).transpose()

#---------plotting-----------
spice_data = orbit_SIII

fig, ax = plt.subplots(2,3)
ax[0,0].plot(spice_data[0], transformed[1], label='calculated', color='r')
ax[0,0].plot(spice_data[0], spice_data[1], label='spice', color='b')
ax[0,0].set_title('X')

ax[0,1].plot(spice_data[0], transformed[2], label='calculated', color='r')
ax[0,1].plot(spice_data[0], spice_data[2], label='spice', color='b')
ax[0,1].set_title('Y')

ax[0,2].plot(spice_data[0], transformed[3], label='calculated', color='r')
ax[0,2].plot(spice_data[0], spice_data[3], label='spice', color='b')
ax[0,2].set_title('Z')

# ax[1,0].plot(transformed[0], transformed[4], label='calculated', color='r')
# ax[1,0].plot(spice_data[0], spice_data[4], label='spice', color='b')
# ax[1,0].set_title('R')

# ax[1,1].plot(transformed[0], transformed[5], label='calculated', color='r')
# ax[1,1].plot(spice_data[0], spice_data[5], label='spice', color='b')
# ax[1,1].set_title('THETA')

# ax[1,2].plot(transformed[0], transformed[6], label='calculated', color='r')
# ax[1,2].plot(spice_data[0], spice_data[6], label='spice', color='b')
# ax[1,2].set_title('PHI')

#ax[0,2].legend()
plt.show()