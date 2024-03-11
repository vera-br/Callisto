# load modules and functions
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d 
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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

time = (orbit_cphio[0] - orbit_CA[0]) / 60

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
phi = 0
real_As = np.linspace(0.25, 1, 4)
phis = np.linspace(0, np.pi/6, 4)
t_longperiod = orbit_cal_JSO_LP[0]
t_LP = (orbit_cal_JSO_LP[0] - orbit_CA[0]) / 60

fig, ax = plt.subplots(3, 1, layout='constrained')

ax[0].plot(time, B_PDS[1], label='Data', color='k', linewidth=0.7, alpha=0.7)
ax[1].plot(time, B_PDS[2], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[2].plot(time, B_PDS[3], label='PDS', color='k', linewidth=0.7, alpha=0.7)

ax[0].plot(time, B_poly[:,0], color='k', linestyle='--', label='Polynomial')
ax[1].plot(time, B_poly[:,1], color='k', linestyle='--')
ax[2].plot(time, B_poly[:,2], color='k', linestyle='--')

colors = ['#ffb000', '#648fff', '#dc267f', '#fe6100']
for real_A, color in zip(real_As, colors):
    aeiphi = real_A * np.exp(-1j * phi)
    B_induced_model_shifted = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, aeiphi=aeiphi, shifted=True, t_longperiod=t_longperiod)
    
    B_total_shifted = B_poly + B_induced_model_shifted

    ax[0].plot(time, B_total_shifted[:, 0], color=color, label='|A| = {}'.format(real_A))
    ax[1].plot(time, B_total_shifted[:, 1], color=color)
    ax[2].plot(time, B_total_shifted[:, 2], color=color)



ax[0].set_ylabel('$B_x$ [nT]')
ax[1].set_ylabel('$B_y$ [nT]')
ax[2].set_ylabel('$B_z$ [nT]')

# ax[0,1].set_ylabel('$B_x$ [nT]')
# ax[1,1].set_ylabel('$B_y$ [nT]')
# ax[2,1].set_ylabel('$B_z$ [nT]')

# ax[0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[2].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

for ax_i in ax.ravel():
    ax_i.tick_params(axis='both', direction='in', top=True, right=True, which='both')
    ax_i.yaxis.set_minor_locator(AutoMinorLocator())

for ax_i  in ax:
    ax_i.set_xlim(min(time), max(time))


ax[0].legend(framealpha=1, fancybox=True, loc='upper right')

titles = ['blah', 'C3', 'C9']
fig.suptitle('Flyby {}'.format(titles[flyby_n]))

plt.show()

