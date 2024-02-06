from trajectories.trajectory_analysis import *
from field_functions import *
from khurana1997 import *
from khurana_2 import *
from jupiter_field import *

#---------constants-----------





#---------load data-----------


galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()

callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'G')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'G')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'G')

flyby_n = 2

orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]

orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag_ = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = orbit_cal_SIII_mag_.copy()

B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)


#--------------------

B_external = Bext_Community(orbit_cal_SIII)
B_external_mag = np.sqrt(B_external[:, 0]**2 + B_external[:, 1]**2 + B_external[:, 2]**2)

fig, ax = plt.subplots(2, 2, figsize=(10, 6))

ax[0, 0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
#ax[0, 0].plot(orbit_cal_SIII[0], B_external[:, 0], label='B_external', color='k')

ax[0, 0].set_title('Bx')
#ax[0, 0].set_xlim(min(orbit_cal_SIII[0]), max(orbit_cal_SIII[0]))

ax[0, 1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
#ax[0, 1].plot(orbit_cal_SIII[0], B_external[:, 1], label='B_external', color='k')

ax[0, 1].set_title('By')
#ax[0, 1].set_xlim(min(orbit_cal_SIII[0]), max(orbit_cal_SIII[0]))

ax[1, 0].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
#ax[1, 0].plot(orbit_cal_SIII[0], B_external[:, 2], label='B_external', color='k')

ax[1, 0].set_title('Bz')
#ax[1, 0].set_xlim(min(orbit_cal_SIII[0]), max(orbit_cal_SIII[0]))

ax[1, 1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)
#ax[1, 1].plot(orbit_cal_SIII[0], B_external_mag, label='B_external', color='k')

ax[1, 1].set_title('|B|')


models = ["Pioneer10", "Voyager1", "Voyager2", "common"]
colour = ["r", "b", "g", "gold"]

for i in range(4):

    B_sheet = B_sheet_khurana2(orbit_cal_JSO, orbit_cal_SIII_mag, models[i])
    B_sheet_mag = np.sqrt(B_sheet[:, 0]**2 + B_sheet[:, 1]**2 + B_sheet[:, 2]**2)

    B_full_ext = B_sheet + B_external
    Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)


    #ax[0, 0].plot(orbit_cal_SIII[0], B_sheet[:, 0], label='B_sheet', color=colour[i], linestyle="--")
    ax[0, 0].plot(orbit_cal_SIII[0], B_full_ext[:, 0], label=models[i], color=colour[i])

    #ax[0, 1].plot(orbit_cal_SIII[0], B_sheet[:, 1], label='B_sheet', color=colour[i], linestyle="--")
    ax[0, 1].plot(orbit_cal_SIII[0], B_full_ext[:, 1], label=models[i], color=colour[i])

    #ax[1, 0].plot(orbit_cal_SIII[0], B_sheet[:, 2], label='B_sheet', color=colour[i], linestyle="--")
    ax[1, 0].plot(orbit_cal_SIII[0], B_full_ext[:, 2], label=models[i], color=colour[i])

    #ax[1, 1].plot(orbit_cal_SIII[0], B_sheet_mag, label='B_sheet', color=colour[i], linestyle="--")
    ax[1, 1].plot(orbit_cal_SIII[0], Bmag_full_ext, label=models[i], color=colour[i])

ax[1, 1].legend()
#ax[1, 1].set_xlim(min(orbit_cal_SIII[0]), max(orbit_cal_SIII[0]))

#plt.suptitle("Flyby C3")
plt.tight_layout()
plt.show()