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
flyby_n = 1

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

J2000 = datetime(2000,1,1,12) # difference between J2000 and UTC


# time = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]



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




B_sheet = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag, orbit_cal_SIII)
B_external = Bext_Community(orbit_SIII)

B_full_ext = B_external + B_sheet
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

B_sheet_cal = B_sheet_khurana(orbit_cal_JSO_LP, orbit_cal_SIII_mag_LP, orbit_cal_SIII_LP)
B_external_cal = Bext_Community(orbit_cal_SIII_LP)

B_full_ext_cal = B_external_cal + B_sheet_cal
t_longperiod = orbit_cal_SIII_LP[0]

conductivities, radii = [], []

real_A = 0.85
phis = np.linspace(0, np.pi/6, 4)
t_longperiod = orbit_cal_JSO_LP[0]

fig, ax = plt.subplots(3, 1, layout='constrained')

colour = 'k'
ax[0].plot(t_longperiod, B_full_ext_cal[:,0], color='k', linestyle=':', label='Full Ext. Cal.')
ax[1].plot(t_longperiod, B_full_ext_cal[:,1], color='k', linestyle=':')
ax[2].plot(t_longperiod, B_full_ext_cal[:,2], color='k', linestyle=':')
ax[0].plot(time, B_PDS[1], label='Data', color='k', linewidth=0.7, alpha=0.7)
ax[1].plot(time, B_PDS[2], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[2].plot(time, B_PDS[3], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[0].plot(time, B_poly[:,0], color=colour, linestyle='--', label='Full Ext.')
ax[1].plot(time, B_poly[:,1], color=colour, linestyle='--')
ax[2].plot(time, B_poly[:,2], color=colour, linestyle='--')

colors = ['k', 'r', 'b', 'g']
for phi, color in zip(phis, colors):
    aeiphi = real_A * np.exp(-1j * phi)
    B_induced_model_shifted = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, aeiphi=aeiphi, shifted=True, t_longperiod=t_longperiod)
    
    B_total_shifted = B_poly + B_induced_model_shifted

    ax[0].plot(time, B_total_shifted[:, 0], color=color, label='$\phi = {}\xb0$'.format(int(np.ceil(phi * 180 / np.pi))))
    ax[1].plot(time, B_total_shifted[:, 1], color=color)
    ax[2].plot(time, B_total_shifted[:, 2], color=color)



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

ax[0].legend(framealpha=1, fancybox=True, loc='upper right')

titles = ['blah', 'C3', 'C9']
fig.suptitle('Flyby {}'.format(titles[flyby_n]))
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

