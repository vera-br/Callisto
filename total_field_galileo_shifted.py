# load modules and functions
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from field_functions import *
from khurana1997 import *
from jupiter_field import *
#from constants import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
callisto_wrt_jupiter_cphio_CA = get_closest_approach_data("callisto", "jupiter", "SIII", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
galileo_wrt_jupiter_SIII_mag = Galileo_trajectories_SIII_mag_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')
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

orbit_cal_SIII_LP = callisto_jupiter_SIII_longperiod["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag_LP = callisto_jupiter_SIII_mag_longperiod["orbit%s" % (flyby_n)]
orbit_cal_JSO_LP = callisto_jupiter_JSO_longperiod["orbit%s" % (flyby_n)]

B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

time = orbit_cphio[0]


#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

#----------polynomial fits-------------
polys = []
section = int(len(orbit_cphio[0]) / 3)
for i in range(3):
    t = orbit_cphio[0]
    t = np.append(t[:section], t[2*section:])
    Bi = B_PDS[i+1]
    Bi = np.append(Bi[:section], Bi[2*section:])
    poly = np.polyfit(t, Bi, 3)
    p = np.poly1d(poly)
    polys.append(p)

B_poly = []
for p in polys:
    Bi_poly = p(orbit_cphio[0])
    B_poly.append(Bi_poly)
B_poly = np.transpose(B_poly)
B_poly_mag = np.sqrt(B_poly[:, 0]**2 + B_poly[:, 1]**2 + B_poly[:, 2]**2)


# induced field parameters
model = 'ocean and iono'

if model == 'ocean and iono':
    # Conducting Ocean and Ionosphere
    r_core = 0.7 * R_C ;   r_ocean = 0.95 * R_C ;   r_surface = R_C    ;   r_iono = 1.042 * R_C
    sig_core = 1e-9    ;   sig_ocean = 3   ;   sig_surface = 1e-9 ;   sig_iono = 0.5e-3

    radii = [r_core, r_ocean, r_surface, r_iono]
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]

elif model == 'ocean only':
    r_core = 0.7 * R_C ;   r_ocean = 0.999 * R_C ;   r_surface = R_C
    sig_core = 1e-9    ;   sig_ocean = 1e0    ;   sig_surface = 1e-9

    radii = [r_core, r_ocean, r_surface]
    conductivities = [sig_core, sig_ocean, sig_surface]

elif model == 'surface ocean':
    r_core = 0.7 * R_C ;   r_ocean = R_C
    sig_core = 1e-9    ;   sig_ocean = 1e0

    radii = [r_core, r_ocean]
    conductivities = [sig_core, sig_ocean]

elif model == 'iono only':
    r_surface = R_C    ;   r_iono = 1.042 * R_C
    sig_surface = 1e-9 ;   sig_iono = 0.5e-3

    radii = [r_surface, r_iono]
    conductivities = [sig_surface, sig_iono]



B_sheet = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag, orbit_cal_SIII)
B_external = Bext_Community(orbit_SIII)

B_full_ext = B_external + B_sheet
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

B_sheet_cal = B_sheet_khurana(orbit_cal_JSO_LP, orbit_cal_SIII_mag_LP, orbit_cal_SIII_LP)
B_external_cal = Bext_Community(orbit_cal_SIII_LP)

B_full_ext_cal = B_external_cal + B_sheet_cal
t_longperiod = orbit_cal_SIII_LP[0]

B_induced_poly = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_full_ext_cal, t_longperiod,  2*np.pi /(10.1*3600), conductivities, radii, Styczinski=True)
Bmag_induced_model_og = np.sqrt(B_induced_poly[:, 0]**2 + B_induced_poly[:, 1]**2 + B_induced_poly[:, 2]**2)
B_total_poly = B_poly + B_induced_poly
B_mag_tot_poly = np.sqrt(B_total_poly[:, 0]**2 + B_total_poly[:, 1]**2 + B_total_poly[:, 2]**2)

#---------plot-----------
# plt.style.use('dark_background')

# fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8.3, 7.7), dpi=300, constrained_layout=True)
fig, ax = plt.subplots(3, 1, sharex=True, constrained_layout=True)

# ax[0].plot(time, Bx_smooth, label='PDS', color='k')
# ax[1].plot(time, By_smooth, label='PDS', color='k')
# ax[2].plot(time, Bz_smooth, label='PDS', color='k')

# ax[0].plot(time, B_PDS[1], label='PDS', color='k', alpha=0.3)
# ax[1].plot(time, B_PDS[2], label='PDS', color='k', alpha=0.3)
# ax[2].plot(time, B_PDS[3], label='PDS', color='k', alpha=0.3)

ax[0].plot(time, B_PDS[1], label='Data', color='k', linewidth=0.7, alpha=0.7)
ax[1].plot(time, B_PDS[2], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[2].plot(time, B_PDS[3], label='PDS', color='k', linewidth=0.7, alpha=0.7)

colour = '#785ef0'
ax[0].plot(time, B_poly[:,0], color=colour, linestyle='--', label='No Induction Model')
ax[1].plot(time, B_poly[:,1], color=colour, linestyle='--')
ax[2].plot(time, B_poly[:,2], color=colour, linestyle='--')

#-----------

colour3 = '#785ef0'
ax[0].plot(time, B_full_ext[:,0], color=colour3, linestyle='-', label='OG. Khur.')
ax[1].plot(time, B_full_ext[:,1], color=colour3, linestyle='-')
ax[2].plot(time, B_full_ext[:,2], color=colour3, linestyle='-')
#----------

ax[0].plot(time, B_total_poly[:, 0], color=colour, label = 'Induction Model')
ax[1].plot(time, B_total_poly[:, 1], color=colour)
ax[2].plot(time, B_total_poly[:, 2], color=colour)

ax[0].set_ylabel('$B_x$ [nT]', fontsize=16)
ax[1].set_ylabel('$B_y$ [nT]', fontsize=16)
ax[2].set_ylabel('$B_z$ [nT]', fontsize=16)

# ax[0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[2].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

ax[0].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
ax[1].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
ax[2].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)

ax[0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
ax[2].yaxis.set_minor_locator(AutoMinorLocator())

ax[0].set_xlim(min(time), max(time))
ax[1].set_xlim(min(time), max(time))
ax[2].set_xlim(min(time), max(time))

# ax[0].legend(framealpha=1, fancybox=True, fontsize=14)

fig.suptitle('Flyby C9', fontsize=16)
# fig.tight_layout()
# plt.savefig('galileo.png', facecolor=('b', 0))
plt.show()


# fig, ax = plt.subplots(2, 2)
# ax[0,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[0,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[1,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[1,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

# ax[0,0].plot(time, Bx_smooth, label='PDS', color='k')
# ax[0,1].plot(time, By_smooth, label='PDS', color='k')
# ax[1,0].plot(time, Bz_smooth, label='PDS', color='k')
# ax[1,1].plot(time, Bmag_smooth, label='PDS Smoothed', color='k')

# ax[0,0].plot(time, B_PDS[1], label='PDS', color='k', alpha=0.3)
# ax[0,1].plot(time, B_PDS[2], label='PDS', color='k', alpha=0.3)
# ax[1,0].plot(time, B_PDS[3], label='PDS', color='k', alpha=0.3)
# ax[1,1].plot(time, B_mag, label='PDS', color='k', alpha=0.3)

# # Polynomial
# ax[0,0].plot(time, B_poly[:,0], '--r')
# ax[0,1].plot(time, B_poly[:,1], '--r')
# ax[1,0].plot(time, B_poly[:,2], '--r')
# ax[1,1].plot(time, B_poly_mag, '--r')

# ax[0,0].plot(time, B_total_poly[:, 0], color='r')
# ax[0,1].plot(time, B_total_poly[:, 1], color='r')
# ax[1,0].plot(time, B_total_poly[:, 2], color='r')
# ax[1,1].plot(time, B_mag_tot_poly, label='Poly.', color='r')

# ax[0,0].set_title('$B_x$')
# ax[0,1].set_title('$B_y$')
# ax[1,0].set_title('$B_z$')
# ax[1,1].set_title('$|B|$')

# ax[0,0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[0,1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[1,0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[1,1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

# ax[0,0].set_xlim(min(time), max(time))
# ax[0,1].set_xlim(min(time), max(time))
# ax[1,0].set_xlim(min(time), max(time))
# ax[1,1].set_xlim(min(time), max(time))

# ax[1,1].legend(framealpha=1, fancybox=True)

# if model == 'ocean and iono':
#     fig.suptitle('$r_{} = {:.1f}-{:.1f} R_C,  r_{} = 1-{:.3f} R_C$ \n $\sigma_{} = {:2.0e} Sm^{}, \sigma_{} = {:2.0e} Sm^{}$'.format('{ocean}', r_core / R_C, r_ocean / R_C, '{iono}', r_iono / R_C, '{ocean}', sig_ocean, '{-1}', '{iono}', sig_iono, '{-1}'))
# elif model == 'ocean only':
#     fig.suptitle('$r_{} = {:.1f}-{:.1f} R_C, \sigma_{} = {:2.0e} Sm^{}$'.format('{ocean}', r_core / R_C, r_ocean / R_C, '{ocean}', sig_ocean, '{-1}'))
# elif model == 'iono only':
#     fig.suptitle('$r_{} = 1-{:.3f} R_C, \sigma_{} = {:2.0e} Sm^{}$'.format('{iono}', r_iono / R_C, '{iono}', sig_iono, '{-1}'))
# plt.show()

