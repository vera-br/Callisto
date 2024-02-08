# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d 
from iminuit import Minuit
from iminuit.cost import LeastSquares

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
callisto_wrt_jupiter_cphio_CA = get_closest_approach_data("callisto", "jupiter", "SIII", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
galileo_wrt_jupiter_SIII_mag = Galileo_trajectories_SIII_mag_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'IAU_JUPITER')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'MAG_VIP4')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = galileo_wrt_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII_CA = callisto_wrt_jupiter_cphio_CA["CA_orbit%s" % (flyby_n)]
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

# SIII_mag coords
B_external = Bext_Community2(orbit_SIII, orbit_SIII_mag)
Bmag_external = np.sqrt(B_external[:, 0]**2 + B_external[:, 1]**2 + B_external[:, 2]**2)
B_sheet = B_sheet_khurana2(orbit_cal_JSO, orbit_SIII_mag)
Bmag_sheet = np.sqrt(B_sheet[:, 0]**2 + B_sheet[:, 1]**2 + B_sheet[:, 2]**2)

B_external_cal = Bext_Community2(orbit_cal_SIII, orbit_cal_SIII_mag)
B_sheet_cal = B_sheet_khurana2(orbit_cal_JSO, orbit_cal_SIII_mag)
B_full_ext_cal = B_external_cal + B_sheet_cal

B_full_ext = B_external + B_sheet
# B_full_ext = convert_B_to_PDS_CPhiO(B_full_ext, orbit_cal_SIII, orbit_cal_SIII_CA)
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

# SIII coords
B_external_og = Bext_Community(orbit_SIII)
Bmag_external_og = np.sqrt(B_external_og[:, 0]**2 + B_external_og[:, 1]**2 + B_external_og[:, 2]**2)
B_sheet_og = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag, orbit_SIII)
Bmag_sheet_og = np.sqrt(B_sheet_og[:, 0]**2 + B_sheet_og[:, 1]**2 + B_sheet_og[:, 2]**2)

B_external_cal_og = Bext_Community(orbit_cal_SIII)
B_sheet_cal_og = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
B_full_ext_cal_og = B_external_cal_og + B_sheet_cal_og

B_full_ext_og = B_external_og + B_sheet_og
# B_full_ext_og = convert_B_to_PDS_CPhiO(B_full_ext_og, orbit_cal_SIII, orbit_cal_SIII_CA)
Bmag_full_ext_og = np.sqrt(B_full_ext_og[:, 0]**2 + B_full_ext_og[:, 1]**2 + B_full_ext_og[:, 2]**2)


#----------polynomial fits-------------
polys = []
section = int(len(B_PDS[0]) / 3)
for i in range(3):
    t = B_PDS[0]
    t = np.append(t[:section], t[2*section:])
    Bi = B_PDS[i+1]
    Bi = np.append(Bi[:section], Bi[2*section:])
    poly = np.polyfit(t, Bi, 3)
    p = np.poly1d(poly)
    polys.append(p)

B_poly = []
for p in polys:
    Bi_poly = p(B_PDS[0])
    B_poly.append(Bi_poly)
B_poly = np.transpose(B_poly)
B_poly_mag = np.sqrt(B_poly[:, 0]**2 + B_poly[:, 1]**2 + B_poly[:, 2]**2)


# note: do minimisation on smoothed PDS minus polyfit background field
def B_minimisation_func_x(t, r_core, r_ocean, r_iono, sig_ocean, sig_iono):
    conductivities = [1e-6, sig_ocean, 1e-6, sig_iono]
    radii = np.array([r_core, r_ocean, 1, r_iono]) * R_C
    p = polys[0]
    B_poly = p(t)
    B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, 2*np.pi /(10.1*3600), conductivities, radii)
    B_total = B_poly + B_induced[:,0]
    return B_total

data_t = B_PDS[0]
data_B = Bx_smooth
data_B_err = 1 * np.ones_like(data_B)
cost = LeastSquares(data_t, data_B, data_B_err, B_minimisation_func_x)
m = Minuit(cost, 0.1, 0.8, 1.1, 1, 0.1)
m.migrad()

# induced field parameters
r_core = 0.1 * R_C ; r_ocean = 0.8 * R_C ; r_surface = R_C    ; r_iono = 1.02 * R_C
sig_core = 1e-6    ; sig_ocean = 1      ; sig_surface = 1e-6 ; sig_iono = 1e-6

radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]


# SIII_mag
# B_induced_model = B_induced_finite_conductivity_multilayer(orbit_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii)
# Bmag_induced_model = np.sqrt(B_induced_model[:, 0]**2 + B_induced_model[:, 1]**2 + B_induced_model[:, 2]**2)
# B_total = B_full_ext + B_induced_model
# B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

# SIII
B_induced_model_og = B_induced_finite_conductivity_multilayer(orbit_cphio, B_full_ext_cal_og, 2*np.pi /(10.1*3600), conductivities, radii)
Bmag_induced_model_og = np.sqrt(B_induced_model_og[:, 0]**2 + B_induced_model_og[:, 1]**2 + B_induced_model_og[:, 2]**2)
B_total_og = B_full_ext_og + B_induced_model_og
B_mag_tot_og = np.sqrt(B_total_og[:, 0]**2 + B_total_og[:, 1]**2 + B_total_og[:, 2]**2)

# Polynomial
B_induced_poly = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, 2*np.pi /(10.1*3600), conductivities, radii)
Bmag_induced_model_og = np.sqrt(B_induced_poly[:, 0]**2 + B_induced_poly[:, 1]**2 + B_induced_poly[:, 2]**2)
B_total_poly = B_poly + B_induced_poly
B_mag_tot_poly = np.sqrt(B_total_poly[:, 0]**2 + B_total_poly[:, 1]**2 + B_total_poly[:, 2]**2)


#---------plot-----------
# plot_time_evolution_Gal(B_total, orbit_cphio, orbit_CA, flyby_n, "Total")

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS Smoothed', color='k')

ax[0,0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax[0,1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
ax[1,0].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
ax[1,1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)

ax[0,0].plot(B_PDS[0], B_poly[:,0], '--k')
ax[0,1].plot(B_PDS[0], B_poly[:,1], '--k')
ax[1,0].plot(B_PDS[0], B_poly[:,2], '--k')
ax[1,1].plot(B_PDS[0], B_poly_mag, '--k')

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

# SIII_mag
# ax[0,0].plot(B_PDS[0], B_full_ext[:,0], '--g', label='Full Ext.')
# ax[0,1].plot(B_PDS[0], B_full_ext[:,1], '--g', label='Full Ext.')
# ax[1,0].plot(B_PDS[0], B_full_ext[:,2], '--g', label='Full Ext.')
# ax[1,1].plot(B_PDS[0], Bmag_full_ext, '--g')

# ax[0,0].plot(orbit_cphio[0], B_total[:, 0], label='Calc.', color='g')
# ax[0,1].plot(orbit_cphio[0], B_total[:, 1], label='Calc.', color='g')
# ax[1,0].plot(orbit_cphio[0], B_total[:, 2], label='Calc.', color='g')
# ax[1,1].plot(orbit_cphio[0], B_mag_tot, label='Calc.', color='g')

# SIII
ax[0,0].plot(B_PDS[0], B_full_ext_og[:,0], '--b', label='Full Ext. OG')
ax[0,1].plot(B_PDS[0], B_full_ext_og[:,1], '--b', label='Full Ext. OG')
ax[1,0].plot(B_PDS[0], B_full_ext_og[:,2], '--b', label='Full Ext. OG')
ax[1,1].plot(B_PDS[0], Bmag_full_ext_og, '--b')

ax[0,0].plot(orbit_cphio[0], B_total_og[:, 0], label='Calc.', color='b')
ax[0,1].plot(orbit_cphio[0], B_total_og[:, 1], label='Calc.', color='b')
ax[1,0].plot(orbit_cphio[0], B_total_og[:, 2], label='Calc.', color='b')
ax[1,1].plot(orbit_cphio[0], B_mag_tot_og, label='Model', color='b')

# Polynomial
ax[0,0].plot(B_PDS[0], B_poly[:,0], '--r', label='Full Ext. OG')
ax[0,1].plot(B_PDS[0], B_poly[:,1], '--r', label='Full Ext. OG')
ax[1,0].plot(B_PDS[0], B_poly[:,2], '--r', label='Full Ext. OG')
ax[1,1].plot(B_PDS[0], B_poly_mag, '--r')

ax[0,0].plot(orbit_cphio[0], B_total_poly[:, 0], label='Calc.', color='r')
ax[0,1].plot(orbit_cphio[0], B_total_poly[:, 1], label='Calc.', color='r')
ax[1,0].plot(orbit_cphio[0], B_total_poly[:, 2], label='Calc.', color='r')
ax[1,1].plot(orbit_cphio[0], B_mag_tot_poly, label='Poly.', color='r')

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

fig, ax = plt.subplots(1,1)
ax.plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax.plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax.plot(B_PDS[0], B_poly[:,0], '--r')
ax.plot(orbit_cphio[0], B_total_poly[:, 0], label='Calc.', color='r')
ax.plot(data_t, B_minimisation_func_x(data_t, *m.values), 'b', label='Minimised')
ax.legend()
plt.show()
