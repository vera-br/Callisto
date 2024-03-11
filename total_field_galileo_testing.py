# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d 
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.optimize import LinearConstraint

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
callisto_jupiter_SIII_longperiod = get_spice_data('callisto', 'jupiter', 'SIII_fullcycle', 'G')
callisto_jupiter_SIII_mag_longperiod = get_spice_data('callisto', 'jupiter', 'SIII_mag_fullcycle', 'G')
callisto_jupiter_JSO_longperiod = get_spice_data('callisto', 'jupiter', 'jupsunorb_fullcycle', 'G')


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

orbit_cal_SIII_LP = callisto_jupiter_SIII_longperiod["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag_LP = callisto_jupiter_SIII_mag_longperiod["orbit%s" % (flyby_n)]
orbit_cal_JSO_LP = callisto_jupiter_JSO_longperiod["orbit%s" % (flyby_n)]
t_longperiod = orbit_cal_SIII_LP[0]

#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

#---------magnetic fields-----------

# SIII coords
B_external_cal = Bext_Community(orbit_cal_SIII_LP)
B_sheet_cal = B_sheet_khurana(orbit_cal_JSO_LP, orbit_cal_SIII_mag_LP, orbit_cal_SIII_LP)
B_full_ext_cal = B_external_cal + B_sheet_cal



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
    conductivities = [1e-9, sig_ocean, 1e-9, sig_iono]
    radii = np.array([r_core, r_ocean, 1, r_iono]) * R_C
    B_induced = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, Styczinski='many', shifted=True, t_longperiod=t_longperiod)
    B_total = B_induced[:,0]
    return B_total

data_t = B_PDS[0]
data_B = Bx_smooth - B_poly[:,0]
data_B_err = 1 * np.ones_like(data_B)
cost = LeastSquares(data_t, data_B, data_B_err, B_minimisation_func_x)

cons = LinearConstraint([[-1,1,0,0,0],[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], [0,0.5,0,1,1e-9,1e-9], [1, 1, 1, 1.1,10,1])
m = Minuit(cost, 0.8, 0.9, 1.1, 1, 0.1)
m.scipy(constraints=cons)
m.migrad()
print(m.values)
# induced field parameters
r_core = m.values['r_core'] * R_C ; r_ocean = m.values['r_ocean'] * R_C ; r_surface = R_C    ; r_iono = m.values['r_iono'] * R_C
sig_core = 1e-9    ; sig_ocean = m.values['sig_ocean']      ; sig_surface = 1e-9 ; sig_iono = m.values['sig_iono']

radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]


# # Polynomial
B_induced_poly = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, Styczinski=True, shifted=True, t_longperiod=t_longperiod)
Bmag_induced_model_og = np.sqrt(B_induced_poly[:, 0]**2 + B_induced_poly[:, 1]**2 + B_induced_poly[:, 2]**2)
B_total_poly = B_poly + B_induced_poly
B_mag_tot_poly = np.sqrt(B_total_poly[:, 0]**2 + B_total_poly[:, 1]**2 + B_total_poly[:, 2]**2)


# #---------plot-----------

fig, ax = plt.subplots(3,1)
ax[0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[2].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
# ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS Smoothed', color='k')

ax[0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax[1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
ax[2].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
# ax[1,1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)

ax[0].plot(B_PDS[0], B_poly[:,0], '--k')
ax[1].plot(B_PDS[0], B_poly[:,1], '--k')
ax[2].plot(B_PDS[0], B_poly[:,2], '--k')
# ax[1,1].plot(B_PDS[0], B_poly_mag, '--k')

ax[0].plot(orbit_cphio[0], B_total_poly[:, 0], label='Calc.', color='r')
ax[1].plot(orbit_cphio[0], B_total_poly[:, 1], label='Calc.', color='r')
ax[2].plot(orbit_cphio[0], B_total_poly[:, 2], label='Calc.', color='r')
# ax[1,1].plot(orbit_cphio[0], B_mag_tot_poly, label='Poly.', color='r')

ax[0].set_title('Bx')
ax[1].set_title('By')
ax[2].set_title('Bz')
# ax[1,1].set_title('|B|')

ax[0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[2].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
# ax[1,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

# ax[1,1].legend(framealpha=1, fancybox=True)

# fig.suptitle('r_ocean = ' + str(r_core / R_C) + '-' + str(r_ocean / R_C) + ' R_C, r_iono = 1-' + str(r_iono / R_C) + ' R_C \n  sig_ocean = ' + str(sig_ocean) + ', sig_iono = ' + str(sig_iono))
plt.show()

# fig, ax = plt.subplots(1,1)
# ax.plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
# ax.plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
# ax.plot(B_PDS[0], B_poly[:,0], '--r')
# ax.plot(orbit_cphio[0], B_total_poly[:, 0], label='Calc.', color='r')
# ax.plot(data_t, B_minimisation_func_x(data_t, *m.values), 'b', label='Minimised')
# ax.legend()
# plt.show()
