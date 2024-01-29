from trajectories.trajectory_analysis import *
from field_functions import *
from khurana1997 import *

callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII1965')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')

flyby_n = 2
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]

orbit_cal_SIII_mag[1:5] = orbit_cal_SIII_mag[1:5] / R_J

orbit_cal_SIII[1:5] = orbit_cal_SIII[1:5] / R_J

orbit_cal_SIII_calc = convert_SIII_mag_to_SIII(orbit_cal_SIII_mag)
# orbit_cal_SIII_mag_calc = convert_SIII_to_SIII_mag(orbit_cal_SIII)


fig, ax = plt.subplots(2,3)
ax[0,0].plot(orbit_cal_SIII[0], orbit_cal_SIII[1])
ax[0,0].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[1])
ax[0,1].plot(orbit_cal_SIII[0], orbit_cal_SIII[2])
ax[0,1].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[2])
ax[0,2].plot(orbit_cal_SIII[0], orbit_cal_SIII[3])
ax[0,2].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[3])
ax[1,0].plot(orbit_cal_SIII[0], orbit_cal_SIII[4])
ax[1,0].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[4])
ax[1,1].plot(orbit_cal_SIII[0], orbit_cal_SIII[5])
ax[1,1].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[5])
ax[1,2].plot(orbit_cal_SIII[0], orbit_cal_SIII[6])
ax[1,2].plot(orbit_cal_SIII_calc[0], orbit_cal_SIII_calc[6])
plt.show()

# fig, ax = plt.subplots(2,3)
# ax[0,0].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[1])
# ax[0,0].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[1])
# ax[0,1].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[2])
# ax[0,1].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[2])
# ax[0,2].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[3])
# ax[0,2].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[3])
# ax[1,0].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[4])
# ax[1,0].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[4])
# ax[1,1].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[5])
# ax[1,1].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[5])
# ax[1,2].plot(orbit_cal_SIII_mag[0], orbit_cal_SIII_mag[6])
# ax[1,2].plot(orbit_cal_SIII_mag_calc[0], orbit_cal_SIII_mag_calc[6])
# plt.show()
