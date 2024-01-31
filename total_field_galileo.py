# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *
from khurana1997 import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII1965')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = convert_orbit_SIII_to_SIII_mag(orbit_SIII)
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

#---------magnetic fields-----------

# jovian field
B_external = Bext_Community2(orbit_cal_SIII, orbit_cal_SIII_mag)
B_external_og = Bext_Community(orbit_cal_SIII)
Bmag_external = np.sqrt(B_external[:, 0]**2 + B_external[:, 1]**2 + B_external[:, 2]**2)
Bmag_external_og = np.sqrt(B_external_og[:, 0]**2 + B_external_og[:, 1]**2 + B_external_og[:, 2]**2)
B_external_cal = Bext_Community2(orbit_cal_SIII, orbit_cal_SIII_mag)
B_external_cal_og = Bext_Community(orbit_cal_SIII)


# current sheet
# B_sheet = B_sheet_Community(orbit_cal_SIII, orbit_cal_SIII_mag)
# B_sheet_cal = B_sheet_Community(orbit_cal_SIII, orbit_cal_SIII_mag)
B_sheet_og = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
# B_sheet_cal_og = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
B_sheet = B_sheet_khurana2(orbit_cal_JSO, orbit_cal_SIII_mag)
B_sheet_cal = B_sheet_khurana2(orbit_cal_JSO, orbit_cal_SIII_mag)
Bmag_sheet = np.sqrt(B_sheet[:, 0]**2 + B_sheet[:, 1]**2 + B_sheet[:, 2]**2)
B_full_ext = B_external + B_sheet
B_full_ext_og = B_external_og + B_sheet_og
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)
Bmag_full_ext_og = np.sqrt(B_full_ext_og[:, 0]**2 + B_full_ext_og[:, 1]**2 + B_full_ext_og[:, 2]**2)

# induced field
# C3
r_core = 0.1 * R_C ; r_ocean = 0.8 * R_C ; r_surface = R_C    ; r_iono = 1.02 * R_C
sig_core = 1e-6    ; sig_ocean = 1      ; sig_surface = 1e-6 ; sig_iono = 1e-6

# C9
# r_core = 0.15 * R_C ; r_ocean = 0.7 * R_C ; r_surface = R_C    ; r_iono = 1.045 * R_C
# sig_core = 1e-6    ; sig_ocean = 1      ; sig_surface = 1e-6 ; sig_iono = 0.01

radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]

B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external_cal + B_sheet_cal, 2*np.pi /(10.1*3600), conductivities, radii)
Bmag_induced = np.sqrt(B_induced[:, 0]**2 + B_induced[:, 1]**2 + B_induced[:, 2]**2)

# total field
B_total = B_external + B_sheet + B_induced
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

#---------plot-----------
plot_time_evolution_Gal(B_total, orbit_cphio, orbit_CA, flyby_n, "Total")

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS Smoothed', color='k')

ax[0,0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax[0,1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
ax[1,0].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
ax[1,1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)

# ax[0,0].plot(B_PDS[0], B_external[:,0], label='Jupiter', color='g')
# ax[0,1].plot(B_PDS[0], B_external[:,1], label='Jupiter', color='g')
# ax[1,0].plot(B_PDS[0], B_external[:,2], label='Jupiter', color='g')
# ax[1,1].plot(B_PDS[0], Bmag_external, label='Jupiter', color='g')

# ax[0,0].plot(B_PDS[0], B_sheet[:,0], label='Sheet', color='m')
# ax[0,1].plot(B_PDS[0], B_sheet[:,1], label='Sheet', color='m')
# ax[1,0].plot(B_PDS[0], B_sheet[:,2], label='Sheet', color='m')
# ax[1,1].plot(B_PDS[0], Bmag_sheet, label='Sheet', color='m')

# ax[0,0].plot(B_PDS[0], B_induced[:,0], label='Induced')
# ax[0,1].plot(B_PDS[0], B_induced[:,1], label='Induced')
# ax[1,0].plot(B_PDS[0], B_induced[:,2], label='Induced')
# ax[1,1].plot(B_PDS[0], Bmag_induced, label='Induced')

ax[0,0].plot(B_PDS[0], B_full_ext[:,0], label='Full Ext.', color='b')
ax[0,1].plot(B_PDS[0], B_full_ext[:,1], label='Full Ext.', color='b')
ax[1,0].plot(B_PDS[0], B_full_ext[:,2], label='Full Ext.', color='b')
ax[1,1].plot(B_PDS[0], Bmag_full_ext, label='Full Ext.', color='b')

ax[0,0].plot(B_PDS[0], B_full_ext_og[:,0], label='Full Ext. OG', color='g')
ax[0,1].plot(B_PDS[0], B_full_ext_og[:,1], label='Full Ext. OG', color='g')
ax[1,0].plot(B_PDS[0], B_full_ext_og[:,2], label='Full Ext. OG', color='g')
ax[1,1].plot(B_PDS[0], Bmag_full_ext_og, label='Full Ext. OG', color='g')

ax[0,0].plot(orbit_cphio[0], B_total[:, 0], label='Calc.', color='r')
ax[0,1].plot(orbit_cphio[0], B_total[:, 1], label='Calc.', color='r')
ax[1,0].plot(orbit_cphio[0], B_total[:, 2], label='Calc.', color='r')
ax[1,1].plot(orbit_cphio[0], B_mag_tot, label='Calc.', color='r')

ax[0,0].set_title('Bx')
ax[0,1].set_title('By')
ax[1,0].set_title('Bz')
ax[1,1].set_title('|B|')

ax[0,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[0,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[1,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[1,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,1].legend(framealpha=1, fancybox=True)



fig.suptitle('r_ocean = ' + str(r_core / R_C) + '-' + str(r_ocean / R_C) + ' R_C, r_iono = 1-' + str(r_iono / R_C) + ' R_C \n  sig_ocean = ' + str(sig_ocean) + ', sig_iono = ' + str(sig_iono))
plt.show()

Bx_smooth_offset = Bx_smooth - 0.5 * (max(Bx_smooth) + min(Bx_smooth)) * np.ones_like(Bx_smooth)
Bx_smooth_norm = Bx_smooth_offset / max(Bx_smooth_offset)

By_smooth_offset = By_smooth - 0.5 * (max(By_smooth) + min(By_smooth)) * np.ones_like(By_smooth)
By_smooth_norm = By_smooth_offset / max(By_smooth_offset)

Bz_smooth_offset = Bz_smooth - 0.5 * (max(Bz_smooth) + min(Bz_smooth)) * np.ones_like(Bz_smooth)
Bz_smooth_norm = Bz_smooth_offset / max(Bz_smooth_offset)

Bmag_smooth_offset = Bmag_smooth - 0.5 * (max(Bmag_smooth) + min(Bmag_smooth)) * np.ones_like(Bmag_smooth)
Bmag_smooth_norm = Bmag_smooth_offset / max(Bmag_smooth_offset)



