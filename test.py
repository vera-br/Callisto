from trajectories.trajectory_analysis import *
from induced_field import *
from field_functions import *
from khurana1997 import *
from khurana_2 import *
from jupiter_field import *

from matplotlib.ticker import AutoMinorLocator
from plot_scripts.plot_field_components import *

#---------constants-----------





#---------load data-----------

juice_wrt_callisto = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')

#---------select flybys-----------

flybys = [6, 15, 16, 17]
flyby_n = 6

orbit_cphio = juice_wrt_callisto["orbit%s" % (flyby_n)]
orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]

orbit_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')


# #---------B-field-----------
B_jupiter = Bext_Community(orbit_SIII)
B_sheet = B_sheet_khurana2(orbit_JSO, orbit_SIII_mag)
B_induced = B_induced_infinite(orbit_cphio, B_jupiter + B_sheet, R_C, R_C-80e3)

B_total = B_jupiter + B_sheet + B_induced


plot_time_evolution(B_total, orbit_cphio, orbit_CA, flyby_n, "External")


#--------------------

fig, ax = plt.subplots(1, 3, figsize=(12, 3))

# Defining and applying the common limits
lim = 6

common_xlim = (-lim, lim)  
common_ylim = (-lim, lim) 

ax[0].set_xlim(common_xlim)
ax[0].set_ylim(common_ylim)
ax[1].set_xlim(common_xlim)
ax[1].set_ylim(common_ylim)
ax[2].set_xlim(common_xlim)
ax[2].set_ylim(common_ylim)

ax[0].grid(color='xkcd:dark blue',alpha =0.1)
ax[1].grid(color='xkcd:dark blue',alpha =0.1)
ax[2].grid(color='xkcd:dark blue',alpha =0.1)

# plot jupiter
jup1 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)
jup2 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)
jup3 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)

ax[0].add_patch(jup1)
ax[1].add_patch(jup2)
ax[2].add_patch(jup3)

ax[0].annotate('Callisto', (1.7, -0.25))
ax[1].annotate('Callisto', (1.7, -0.25))
ax[2].annotate('Callisto', (1.7, -0.25))

# set axis labels
ax[0].set_xlabel('x [$R_C$]')
ax[0].set_ylabel('y [$R_C$]')
ax[1].set_xlabel('z [$R_C$]')
ax[1].set_ylabel('y [$R_C$]')
ax[2].set_xlabel('x [$R_C$]')
ax[2].set_ylabel('z [$R_C$]')

ax[0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[2].tick_params(axis='both', direction='in',top = True, right = True, which='both')

ax[0].xaxis.set_minor_locator(AutoMinorLocator()) 
ax[1].xaxis.set_minor_locator(AutoMinorLocator()) 
ax[2].xaxis.set_minor_locator(AutoMinorLocator()) 

ax[0].yaxis.set_minor_locator(AutoMinorLocator()) 
ax[1].yaxis.set_minor_locator(AutoMinorLocator()) 
ax[2].yaxis.set_minor_locator(AutoMinorLocator()) 


# plot closest appproaches
for i in range(len(flybys)):

    flyby_n = flybys[i]
    orbit = juice_wrt_callisto["orbit%s" % (flyby_n)]

    # rho = np.sqrt((orbit_SIII_mag[1])**2 + (orbit_SIII_mag[2]**2))

    ax[0].plot(orbit[1] / R_C, orbit[2] / R_C, label="orbit%s" % (flyby_n))
    ax[2].plot(orbit[1] / R_C, orbit[3] / R_C, label="orbit%s" % (flyby_n))
    ax[1].plot(orbit[3] / R_C, orbit[2] / R_C, label="orbit%s" % (flyby_n))
    

# ax[2].axhspan(2.5, -2.5, facecolor='slateblue', alpha=0.5, label="plasma sheet")

ax[2].legend(loc='center left', bbox_to_anchor=(1.1, 0.5))

plt.tight_layout()
plt.show()