# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
from pandas import Timestamp

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
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

Bx_mins = [-10.75,   2.75,  -10.5,    -2]
Bx_maxs = [  -2.5,   8.25,   -2.5,    15]

By_mins = [   -45,     10,    -54,    21]
By_maxs = [   -20,     40,    -25,    47]

Bz_mins = [   -14,  -22.5,    -19,   -19]
Bz_maxs = [    14,   1.75,     17,  1.75]

flyby_ns = [6, 15, 16, 17]


# specify orbit
i = 0
for flyby_n in flyby_ns:

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
    time_min = Timestamp((time_CA - timedelta(seconds=1800)).strftime('%Y-%m-%d %H:%M:%S'))
    time_max = Timestamp((time_CA + timedelta(seconds=1800)).strftime('%Y-%m-%d %H:%M:%S'))

    ocean_conductivities = [5, 0.01, 0.001]
    colours = ['#785ef0', '#fe6100', '#dc267f']

    #---------inducing field-----------
    B_external = Bext_Community(orbit_juice_SIII)
    B_sheet = B_sheet_khurana(orbit_juice_JSO, orbit_juice_SIII_mag, orbit_cal_SIII)
    B_full_ext = B_external + B_sheet

    B_external_cal = Bext_Community(orbit_cal_SIII)
    B_sheet_cal = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)
    B_full_ext_cal = B_external_cal + B_sheet_cal

    Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8.3, 7.7), dpi=300, constrained_layout=True)

    ax[0].axvline(x=time_CA, color='k', linestyle=":", zorder=0, alpha=0.4)
    ax[1].axvline(x=time_CA, color='k', linestyle=":", zorder=0, alpha=0.4)
    ax[2].axvline(x=time_CA, color='k', linestyle=":", zorder=0, alpha=0.4)

    ax[0].plot(time, B_full_ext[:,0], color='k', linestyle='--', label='Full Ext.', alpha=0.7)
    ax[1].plot(time, B_full_ext[:,1], color='k', linestyle='--', alpha=0.7)
    ax[2].plot(time, B_full_ext[:,2], color='k', linestyle='--', alpha=0.7)

    for sig_ocean_i, colour in zip(ocean_conductivities, colours):
    
        radii_ocean_iono = [0.7*R_C, 0.95*R_C, R_C, 1.042*R_C]
        radii_ocean = [0.7*R_C, 0.95*R_C, R_C]
        sigs_ocean_iono = [1e-9, sig_ocean_i, 1e-9, 0.5e-3]
        sigs_ocean = [1e-9, sig_ocean_i, 1e-9]
    
        # induced field - ocean and iono

        B_induced_model = B_induced_finite_conductivity_multilayer(orbit_juice_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), sigs_ocean_iono, radii_ocean_iono)
        Bmag_induced_model = np.sqrt(B_induced_model[:, 0]**2 + B_induced_model[:, 1]**2 + B_induced_model[:, 2]**2)
        B_total = B_full_ext + B_induced_model
        B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

        ax[0].plot(time, B_total[:, 0], color=colour, linestyle='-', label='$\sigma_{}$ = {}'.format('{ocean}',sig_ocean_i))
        ax[1].plot(time, B_total[:, 1], color=colour, linestyle='-')
        ax[2].plot(time, B_total[:, 2], color=colour, linestyle='-')

        # induced field - ocean and iono

        B_induced_model = B_induced_finite_conductivity_multilayer(orbit_juice_cphio, B_full_ext_cal, 2*np.pi /(10.1*3600), sigs_ocean, radii_ocean)
        Bmag_induced_model = np.sqrt(B_induced_model[:, 0]**2 + B_induced_model[:, 1]**2 + B_induced_model[:, 2]**2)
        B_total = B_full_ext + B_induced_model
        B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

        ax[0].plot(time, B_total[:, 0], color=colour, linestyle='--')
        ax[1].plot(time, B_total[:, 1], color=colour, linestyle='--')
        ax[2].plot(time, B_total[:, 2], color=colour, linestyle='--')

    ax[0].set_ylabel('$B_x$ [nT]', fontsize=16)
    ax[1].set_ylabel('$B_y$ [nT]', fontsize=16)
    ax[2].set_ylabel('$B_z$ [nT]', fontsize=16)

    ax[0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
    ax[1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
    ax[2].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

    ax[0].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
    ax[1].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
    ax[2].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)

    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[2].yaxis.set_minor_locator(AutoMinorLocator())

    ax[0].set_xlim(time_min, time_max)
    ax[1].set_xlim(time_min, time_max)
    ax[2].set_xlim(time_min, time_max)

    ax[0].set_ylim(Bx_mins[i], Bx_maxs[i])
    ax[1].set_ylim(By_mins[i], By_maxs[i])
    ax[2].set_ylim(Bz_mins[i], Bz_maxs[i])

    ax[0].legend(framealpha=1, fancybox=True)

    fig.suptitle('Flyby C{}'.format(flyby_n), fontsize=16)
    i += 1
    plt.savefig('juice{}.png'.format(flyby_n), facecolor=('b',0))
    plt.show()