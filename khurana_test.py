from trajectories.trajectory_analysis import *
from field_functions import *
from khurana1997 import *
from jupiter_field import *
from current_sheet import *
from scipy.ndimage import uniform_filter1d 

galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII1965')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')

flyby_n = 3
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag_gal = convert_orbit_SIII_to_SIII_mag(orbit_SIII)
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

B_sheet = B_sheet_khurana2(orbit_cal_JSO, orbit_SIII_mag_gal)
B_sheet_mag = np.sqrt(B_sheet[:,0]**2 + B_sheet[:,1]**2 + B_sheet[:,2]**2) 

B_sheet_com = B_sheet_Community(orbit_SIII)
B_sheet_com_mag = np.sqrt(B_sheet_com[:,0]**2 + B_sheet_com[:,1]**2 + B_sheet_com[:,2]**2)

B_external = Bext_Community2(orbit_SIII, orbit_SIII_mag_gal)
Bmag_external = np.sqrt(B_external[:,0]**2 + B_external[:,1]**2, B_external[:,2]**2)

B_full_ext = B_sheet + B_external
B_full_ext_com = B_sheet_com + B_external

Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)
Bmag_full_ext_com = np.sqrt(B_full_ext_com[:, 0]**2 + B_full_ext_com[:, 1]**2 + B_full_ext_com[:, 2]**2)

fig, ax = plt.subplots(2,2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,0].plot(B_PDS[0], B_full_ext[:,0], label='Khur.', color='b')
# ax[0,0].plot(B_PDS[0], B_full_ext_com[:,0], label='Comm.', color='r')
# ax[0,0].plot(B_PDS[0], B_sheet[:,0],  '--b', label='Khur.')
# ax[0,0].plot(B_PDS[0], B_sheet_com[:,0],  '--r', label='Comm.')
# ax[0,0].plot(B_PDS[0], B_external[:,0],  '--g', label='Ext.')
ax[0,0].set_title('Bx')
ax[0,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], B_full_ext[:,1], label='Khur.', color='b')
# ax[0,1].plot(B_PDS[0], B_full_ext_com[:,1], label='Comm.', color='r')
# ax[0,1].plot(B_PDS[0], B_sheet[:,1], '--b', label='Khur.')
# ax[0,1].plot(B_PDS[0], B_sheet_com[:,1], '--r', label='Comm.')
# ax[0,1].plot(B_PDS[0], B_external[:,1],  '--g', label='Ext.')
ax[0,1].set_title('By')
ax[0,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], B_full_ext[:,2], label='Khur.', color='b')
# ax[1,0].plot(B_PDS[0], B_full_ext_com[:,2], label='Comm.', color='r')
# ax[1,0].plot(B_PDS[0], B_sheet[:,2], '--b', label='Khur.')
# ax[1,0].plot(B_PDS[0], B_sheet_com[:,2], '--r', label='Comm.')
# ax[1,0].plot(B_PDS[0], B_external[:,2],  '--g', label='Ext.')
ax[1,0].set_title('Bz')
ax[1,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_full_ext, label='Khur.', color='b')
# ax[1,1].plot(B_PDS[0], Bmag_full_ext_com, label='Comm.', color='r')
# ax[1,1].plot(B_PDS[0], B_sheet_mag, '--b', label='Khur.')
# ax[1,1].plot(B_PDS[0], B_sheet_com_mag, '--r', label='Comm.')
# ax[1,1].plot(B_PDS[0], Bmag_external,  '--g', label='Ext.')
ax[1,1].set_title('|B|')
ax[1,1].legend()
ax[1,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

plt.show()
