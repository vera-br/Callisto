from trajectories.trajectory_analysis import *
from field_functions import *
from khurana1997 import *
from jupiter_field import *
from scipy.ndimage import uniform_filter1d 

galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')

flyby_n = 1
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag_gal = convert_SIII_to_SIII_mag(orbit_SIII)
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

B_sheet = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag_gal, orbit_SIII)
B_external = Bext_Community(orbit_cal_SIII)
B_full_ext = B_sheet + B_external
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

fig, ax = plt.subplots(2,2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,0].plot(B_PDS[0], B_full_ext[:,0], label='Full Ext.', color='b')
ax[0,0].set_title('Bx')
ax[0,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], B_full_ext[:,1], label='Full Ext.', color='b')
ax[0,1].set_title('By')
ax[0,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], B_full_ext[:,2], label='Full Ext.', color='b')
ax[1,0].set_title('Bz')
ax[1,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_full_ext, label='Full Ext.', color='b')
ax[1,1].set_title('|B|')
ax[1,1].legend()
ax[1,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

plt.show()
