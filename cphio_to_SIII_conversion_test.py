from trajectories.trajectory_analysis import *
from field_functions import *

callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
juice_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_jupiter_SIII = get_spice_data('juice', 'jupiter', 'SIII', 'J')

flyby_n = 1

O_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
O_cphio = juice_callisto_cphio["orbit%s" % (flyby_n)]
O_SIII = juice_jupiter_SIII["orbit%s" % (flyby_n)]
O_SIII_calc = convert_CPhiO_to_SIII(O_cal_SIII, O_cphio)

O_SIII[1:5] /= R_J
O_SIII_calc[1:5] /= R_J

fig, ax = plt.subplots(2,3)
ax[0,0].plot(O_SIII[0], O_SIII[1])
ax[0,0].plot(O_SIII_calc[0], O_SIII_calc[1])
ax[0,1].plot(O_SIII[0], O_SIII[2])
ax[0,1].plot(O_SIII_calc[0], O_SIII_calc[2])
ax[0,2].plot(O_SIII[0], O_SIII[3])
ax[0,2].plot(O_SIII_calc[0], O_SIII_calc[3])
ax[1,0].plot(O_SIII[0], O_SIII[4])
ax[1,0].plot(O_SIII_calc[0], O_SIII_calc[4])
ax[1,1].plot(O_SIII[0], O_SIII[5])
ax[1,1].plot(O_SIII_calc[0], O_SIII_calc[5])
ax[1,2].plot(O_SIII[0], O_SIII[6])
ax[1,2].plot(O_SIII_calc[0], O_SIII_calc[6])
plt.show()