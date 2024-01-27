from trajectories.trajectory_analysis import *
from field_functions import *
from khurana1997 import *
from khurana_2 import *
from jupiter_field import *
from scipy.ndimage import uniform_filter1d 

galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')

flyby_n = 2
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)


#--------------------

B_Archie = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag)
B_Vera = B_khurana_2(orbit_cal_JSO, orbit_cal_SIII_mag)

fig, ax = plt.subplots(2,2)
ax[0,0].plot(B_PDS[0], B_Archie[0], label='A', color='k')
ax[0,0].plot(B_PDS[0], B_Vera[0], label='V', color='b')
ax[0,0].set_title('Brho')

ax[0,1].plot(B_PDS[0], B_Archie[1], label='A', color='k')
ax[0,1].plot(B_PDS[0], B_Vera[1], label='V', color='b')
ax[0,1].set_title('Bphi')

ax[1,0].plot(B_PDS[0], B_Archie[2], label='A', color='k')
ax[1,0].plot(B_PDS[0], B_Vera[2], label='V', color='b')
ax[1,0].set_title('Bz')

plt.show()