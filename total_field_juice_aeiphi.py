# load modules and functions
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
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
juice_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_jupiter_SIII = get_spice_data('juice', 'jupiter', 'SIII', 'J')
juice_jupiter_SIII_mag = get_spice_data('juice', 'jupiter', 'SIII_mag', 'J')
juice_jupiter_JSO = get_spice_data('juice', 'jupiter', 'jupsunorb', 'J')
juice_callisto_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')


# specify orbit
# for closest approach altitude and SIII_mag z - C9 ~ 17, C3 ~ 19
# 7 = max above, 13 = max below, 11 = closest to 0
# 5 - z = 1.5, 15 - z = -1.5, 6 - z = 2.7, 3 - z = -2.8
flyby_n = 4

orbit_juice_cphio = juice_callisto_cphio["orbit%s" % (flyby_n)]
orbit_juice_SIII = juice_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_juice_SIII_mag = juice_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_juice_JSO = juice_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_juice_CA = juice_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

time = (orbit_juice_cphio[0] - orbit_juice_CA[0]) / 60
mask = (time > -45) & (time < 45)
#---------magnetic fields-----------

# SIII coords
B_external = Bext_Community(orbit_juice_SIII)
B_sheet = B_sheet_khurana(orbit_juice_JSO, orbit_juice_SIII_mag, orbit_cal_SIII)
B_full_ext = B_external + B_sheet

B_external_cal = Bext_Community(orbit_cal_SIII)
B_sheet_cal = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
B_full_ext_cal = B_external_cal + B_sheet_cal

Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

conductivities, radii = [], []

real_A = 0.85
phis = np.linspace(0, np.pi/6, 4)

fig, ax = plt.subplots(3, 2)

colors = ['k', '#648fff', '#dc267f', '#fe6100']
linestyles = ['-', '--', '-.', ':']
for phi, color in zip(phis, colors):
    aeiphi = real_A * np.exp(-1j * phi)
    # B_induced_model = B_induced_finite_conductivity_multilayer(orbit_juice_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, aeiphi=aeiphi)
    B_induced_model_shifted = B_induced_finite_conductivity_multilayer(orbit_juice_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii, aeiphi=aeiphi, shifted=True)
    # Bmag_induced_model = np.sqrt(B_induced_model[:, 0]**2 + B_induced_model[:, 1]**2 + B_induced_model[:, 2]**2)
    # B_total = B_full_ext + B_induced_model
    B_total_shifted = B_full_ext + B_induced_model_shifted
    # B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

    # ax[0].plot(time, B_total[:, 0], color=color, linestyle='--')
    # ax[1].plot(time, B_total[:, 1], color=color, linestyle='--')
    # ax[2].plot(time, B_total[:, 2], linestyle='--')
    ax[0,0].plot(time, B_total_shifted[:, 0], color=color, label='$\phi = {}\xb0$'.format(int(np.ceil(phi * 180 / np.pi))))
    ax[1,0].plot(time, B_total_shifted[:, 1], color=color)
    ax[2,0].plot(time, B_total_shifted[:, 2], color=color)

    phi_t = 10.18 * 60 / 360 * np.ceil(phi * 180 / np.pi)
    ax[0,0].axvline(-phi_t, color=color, linestyle='--', alpha=0.3)
    ax[1,0].axvline(-phi_t, color=color, linestyle='--', alpha=0.3)
    ax[2,0].axvline(-phi_t, color=color, linestyle='--', alpha=0.3)

    ax[0,1].plot(time[mask], B_total_shifted[:, 0][mask], color=color, label='$\phi = {}\xb0$'.format(int(np.ceil(phi * 180 / np.pi))))
    ax[1,1].plot(time[mask], B_total_shifted[:, 1][mask], color=color)
    ax[2,1].plot(time[mask], B_total_shifted[:, 2][mask], color=color)
    

#---------plot-----------
colour = 'k'
ax[0,0].plot(time, B_full_ext_cal[:,0], color='k', linestyle=':')
ax[1,0].plot(time, B_full_ext_cal[:,1], color='k', linestyle=':')
ax[2,0].plot(time, B_full_ext_cal[:,2], color='k', linestyle=':')
ax[0,0].plot(time, B_full_ext[:,0], color=colour, linestyle='--')
ax[1,0].plot(time, B_full_ext[:,1], color=colour, linestyle='--')
ax[2,0].plot(time, B_full_ext[:,2], color=colour, linestyle='--')

# ax[0,1].plot(time[mask], B_full_ext_cal[:,0][mask], color='k', linestyle=':', label='Full Ext. Cal.')
# ax[1,1].plot(time[mask], B_full_ext_cal[:,1][mask], color='k', linestyle=':')
# ax[2,1].plot(time[mask], B_full_ext_cal[:,2][mask], color='k', linestyle=':')
ax[0,1].plot(time[mask], B_full_ext[:,0][mask], color=colour, linestyle='--')
ax[1,1].plot(time[mask], B_full_ext[:,1][mask], color=colour, linestyle='--')
ax[2,1].plot(time[mask], B_full_ext[:,2][mask], color=colour, linestyle='--')


ax[0,0].set_ylabel('$B_x$ [nT]')
ax[1,0].set_ylabel('$B_y$ [nT]')
ax[2,0].set_ylabel('$B_z$ [nT]')


for ax_i in ax.ravel():
    ax_i.tick_params(axis='both', direction='in', top=True, right=True, which='both')
    ax_i.yaxis.set_minor_locator(AutoMinorLocator())
for ax_i in ax[:,0]:
    ax_i.set_xlim(-180, 180)
# print(time[mask])
for ax_i in ax[:,1]:
    ax_i.set_xlim(np.min(time[mask]), np.max(time[mask]))


ax[0,1].legend(framealpha=1, fancybox=True, loc='upper right')

mark_inset(ax[0,0], ax[0,1], loc1=2, loc2=3, linestyle='-', alpha=0.5)
mark_inset(ax[1,0], ax[1,1], loc1=2, loc2=3, linestyle='-', alpha=0.5)
mark_inset(ax[2,0], ax[2,1], loc1=2, loc2=3, linestyle='-', alpha=0.5)

fig.suptitle('Flyby C{}'.format(flyby_n))
plt.show()



