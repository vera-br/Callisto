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
juice_callisto_cphio = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_jupiter_SIII = get_spice_data('juice', 'jupiter', 'SIII', 'J')
juice_jupiter_SIII_mag = get_spice_data('juice', 'jupiter', 'SIII_mag', 'J')
juice_jupiter_JSO = get_spice_data('juice', 'jupiter', 'jupsunorb', 'J')
juice_callisto_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')


# specify orbit
flyby_n = 10

orbit_juice_cphio = juice_callisto_cphio["orbit%s" % (flyby_n)]
orbit_juice_SIII = juice_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_juice_SIII_mag = juice_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_juice_JSO = juice_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_juice_CA = juice_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

J2000 = datetime(2000,1,1,12) # difference between J2000 and UTC
time = orbit_juice_cphio[0]

time = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]
time_CA = Timestamp((J2000 + timedelta(seconds=orbit_juice_CA[0])).strftime('%Y-%m-%d %H:%M:%S'))

#---------magnetic fields-----------

# SIII coords
B_external = Bext_Community(orbit_juice_SIII)
B_sheet = B_sheet_khurana(orbit_juice_JSO, orbit_juice_SIII_mag, orbit_cal_SIII)
B_full_ext = B_external + B_sheet

B_external_cal = Bext_Community(orbit_cal_SIII)
B_sheet_cal = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
B_full_ext_cal = B_external_cal + B_sheet_cal

Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

# induced field parameters
model = 'ocean and iono'

if model == 'ocean and iono':
    # Conducting Ocean and Ionosphere
    r_core = 0.7 * R_C ;   r_ocean = 0.95 * R_C ;   r_surface = R_C    ;   r_iono = 1.042 * R_C
    sig_core = 1e-9    ;   sig_ocean = 5   ;   sig_surface = 1e-9 ;   sig_iono = 0.5e-3

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

B_induced_model = B_induced_finite_conductivity_multilayer(orbit_juice_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), conductivities, radii)
Bmag_induced_model = np.sqrt(B_induced_model[:, 0]**2 + B_induced_model[:, 1]**2 + B_induced_model[:, 2]**2)
B_total = B_full_ext + B_induced_model
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)



#---------plot-----------
plt.style.use('dark_background')

fig, ax = plt.subplots(3, 1)

colour = '#785ef0'

ax[0].plot(time, B_full_ext[:,0], color=colour, linestyle='--', label='Full Ext.')
ax[1].plot(time, B_full_ext[:,1], color=colour, linestyle='--')
ax[2].plot(time, B_full_ext[:,2], color=colour, linestyle='--')

ax[0].plot(time, B_total[:, 0], color=colour, label='Total')
ax[1].plot(time, B_total[:, 1], color=colour)
ax[2].plot(time, B_total[:, 2], color=colour)

ax[0].set_ylabel('$B_x$ [nT]')
ax[1].set_ylabel('$B_y$ [nT]')
ax[2].set_ylabel('$B_z$ [nT]')

ax[0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
ax[1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
ax[2].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

ax[0].tick_params(axis='both', direction='in', top=True, right=True, which='both')
ax[1].tick_params(axis='both', direction='in', top=True, right=True, which='both')
ax[2].tick_params(axis='both', direction='in', top=True, right=True, which='both')

ax[0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
ax[2].yaxis.set_minor_locator(AutoMinorLocator())

ax[0].set_xlim(min(time), max(time))
ax[1].set_xlim(min(time), max(time))
ax[2].set_xlim(min(time), max(time))

ax[0].legend(framealpha=1, fancybox=True)

plt.show()

# fig, ax = plt.subplots(2, 2)

# ax[0,0].plot(time, B_full_ext[:,0], '--g', label='Full Ext.')
# ax[0,1].plot(time, B_full_ext[:,1], '--g', label='Full Ext.')
# ax[1,0].plot(time, B_full_ext[:,2], '--g', label='Full Ext.')
# ax[1,1].plot(time, Bmag_full_ext, '--g')

# ax[0,0].plot(time, B_total[:, 0], label='Calc.', color='g')
# ax[0,1].plot(time, B_total[:, 1], label='Calc.', color='g')
# ax[1,0].plot(time, B_total[:, 2], label='Calc.', color='g')
# ax[1,1].plot(time, B_mag_tot, label='Calc.', color='g')

# ax[0,0].set_title('Bx')
# ax[0,1].set_title('By')
# ax[1,0].set_title('Bz')
# ax[1,1].set_title('|B|')

# ax[0,0].set_xlim(min(time), max(time))
# ax[0,1].set_xlim(min(time), max(time))
# ax[1,0].set_xlim(min(time), max(time))
# ax[1,1].set_xlim(min(time), max(time))

# ax[1,1].legend(framealpha=1, fancybox=True)

# plt.show()

