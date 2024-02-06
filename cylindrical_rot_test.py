from trajectories.trajectory_analysis import *

juice_jupiter_cphio = get_spice_data('juice', 'jupiter', 'cphio', 'J')
juice_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_jupiter_SIII_mag = get_spice_data('juice', 'jupiter', 'SIII_mag', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')
callisto_jupiter_cphio = get_spice_data('callisto', 'jupiter', 'cphio', 'J')

flyby_n = 10

orbit_juice_jupiter_cphio = juice_jupiter_cphio["orbit%s" % (flyby_n)]
orbit_juice_cal_cphio = juice_callisto_cphio["orbit%s" % (flyby_n)]
orbit_juice_jupiter_SIII_mag = juice_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_jupiter_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_jupiter_cphio = callisto_jupiter_cphio["orbit%s" % (flyby_n)]
O_juice_jupiter_SIII_mag = orbit_juice_jupiter_SIII_mag.copy()
O_juice_jupiter_cphio = orbit_juice_jupiter_cphio.copy()
O_cal_jupiter_cphio = orbit_cal_jupiter_cphio.copy()
O_juice_cal_cphio = orbit_juice_cal_cphio.copy()

theta_SIII_mag = orbit_cal_jupiter_SIII_mag[5]
phi_SIII_mag = orbit_cal_jupiter_SIII_mag[6]

y_cal_jup = orbit_cal_jupiter_cphio[2]

O_juice_jupiter_SIII_mag_xyz = O_juice_jupiter_SIII_mag[1:4].transpose()
O_JJ_CPhiO_xyz_calc = []
for phi_i, theta_i, O_JJ_S3M_xyz, y_cal_jup_i in zip(phi_SIII_mag, theta_SIII_mag, O_juice_jupiter_SIII_mag_xyz, y_cal_jup):
    rot_matrix_phi = [[   np.cos(phi_i),   np.sin(phi_i),   0],
                      [  -np.sin(phi_i),   np.cos(phi_i),   0],
                      [               0,               0,   1]]
    rot_matrix_theta = [[ np.sin(theta_i),   0,   np.cos(theta_i)], 
                        [               0,   1,                 0], 
                        [-np.cos(theta_i),   0,   np.sin(theta_i)]]
    O_JJ_CPhiO_calc_i = np.dot(rot_matrix_theta, np.dot(rot_matrix_phi, O_JJ_S3M_xyz))
    x, y, z = O_JJ_CPhiO_calc_i
    O_JJ_CPhiO_calc_i = [y, -x - y_cal_jup_i, z]
    O_JJ_CPhiO_xyz_calc.append(O_JJ_CPhiO_calc_i)

O_JJ_CPhiO_xyz_calc = np.transpose(O_JJ_CPhiO_xyz_calc) / R_C
O_juice_cal_cphio[1:4] = O_juice_cal_cphio[1:4] / R_C

fig, ax = plt.subplots(1,3)
ax[0].plot(O_juice_jupiter_cphio[0], O_juice_cal_cphio[1])
ax[1].plot(O_juice_jupiter_cphio[0], O_juice_cal_cphio[2])
ax[2].plot(O_juice_jupiter_cphio[0], O_juice_cal_cphio[3], label='Spice')

ax[0].plot(O_juice_jupiter_cphio[0], O_JJ_CPhiO_xyz_calc[0])
ax[1].plot(O_juice_jupiter_cphio[0], O_JJ_CPhiO_xyz_calc[1])
ax[2].plot(O_juice_jupiter_cphio[0], O_JJ_CPhiO_xyz_calc[2], label='Calc.')

ax[2].legend()

plt.show()




    
    

