# plotting juice predictions in a plot showing timeseries and CphiO coords

#---------import functions-----------

from trajectories.trajectory_analysis import *
from induced_field import *
from field_functions import *
from khurana1997 import *
from khurana_2 import *
from jupiter_field import *

from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec

#---------load data-----------

juice_wrt_callisto = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')



#---------create figure layout-----------

# Create subplots with 4 rows and 2 columns
fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(11, 6), sharex=True)

# clear figure layout
plt.clf()

# Specify the width ratios of the columns
gs = gridspec.GridSpec(4, 2, width_ratios=[3, 1])

# create subplots for timetraces
ax_l1 = plt.subplot(gs[0, 0])
ax_l2 = plt.subplot(gs[1, 0])
ax_l3 = plt.subplot(gs[2, 0])
ax_l4 = plt.subplot(gs[3, 0])

# Create subplots for the trajectories
ax_r1 = plt.subplot(gs[0:2, 1])
ax_r2 = plt.subplot(gs[2:4, 1])


#---------formatting of timetraces-----------

# add a vertical line at CA
ax_l1.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
ax_l2.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
ax_l3.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
ax_l4.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)

# Apply limits on the time axis
lim = 45
time_lim = (-lim, lim)  

ax_l1.set_xlim(time_lim)
ax_l2.set_xlim(time_lim)
ax_l3.set_xlim(time_lim)
ax_l4.set_xlim(time_lim)

# set axis labels
ax_l1.set_ylabel("$B_x$ [nT]")
ax_l2.set_ylabel("$B_y$ [nT]")
ax_l3.set_ylabel("$B_z$ [nT]")
ax_l4.set_ylabel("|B| [nT]")
ax_l4.set_xlabel("Time after CA [min]")


#---------formatting of trajectories-----------

# plot Callisto
cal1 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)
cal2 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)

ax_r1.add_patch(cal1)
ax_r2.add_patch(cal2)

ax_r1.annotate('Callisto', (1.7, -0.25))
ax_r2.annotate('Callisto', (1.7, -0.25))

# Apply limits on the axis
lim = 6
orbit_lim = (-lim, lim) 

ax_r1.set_xlim(orbit_lim)
ax_r2.set_xlim(orbit_lim)
ax_r1.set_ylim(orbit_lim)
ax_r2.set_ylim(orbit_lim)

# set axis labels
# ax_r1.set_xlabel('x [$R_C$]')
ax_r2.set_xlabel('x [$R_C$]')
ax_r1.set_ylabel('z [$R_C$]')
ax_r2.set_ylabel('y [$R_C$]')


#---------general formatting-----------

# set axis ticks
ax_l1.tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax_l2.tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax_l3.tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax_l4.tick_params(axis='both', direction='in',top = True, right = True, which='both')

ax_r1.tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax_r2.tick_params(axis='both', direction='in',top = True, right = True, which='both')


ax_l1.xaxis.set_minor_locator(AutoMinorLocator())
ax_l2.xaxis.set_minor_locator(AutoMinorLocator())
ax_l3.xaxis.set_minor_locator(AutoMinorLocator())
ax_l4.xaxis.set_minor_locator(AutoMinorLocator())

ax_r1.xaxis.set_minor_locator(AutoMinorLocator()) 
ax_r2.xaxis.set_minor_locator(AutoMinorLocator()) 


ax_l1.yaxis.set_minor_locator(AutoMinorLocator())
ax_l2.yaxis.set_minor_locator(AutoMinorLocator())
ax_l3.yaxis.set_minor_locator(AutoMinorLocator())
ax_l4.yaxis.set_minor_locator(AutoMinorLocator())

ax_r1.yaxis.set_minor_locator(AutoMinorLocator()) 
ax_r2.yaxis.set_minor_locator(AutoMinorLocator()) 



#---------select flybys-----------

flybys = [16, 17, 19] # [3, 4, 6, 7, 10, 13, 16, 17, 19]

surface_layer = 80e3 #m
conductivities=[0.1, 4, 0.1, 1]
radii=[0.6*R_C, R_C-surface_layer, R_C, R_C+100e3]

for flyby_n in flybys:
    
    orbit_cphio = juice_wrt_callisto["orbit%s" % (flyby_n)]
    orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
    orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
    orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]

    orbit_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')
    orbit_CA = orbit_CA["CA_orbit%s" % (flyby_n)]

    # convert time series to minutes before/after CA
    minutes = [(seconds - orbit_CA[0])/60 for seconds in orbit_cphio[0]]

    #---------compute B-field-----------
    B_jupiter = Bext_Community(orbit_SIII)
    B_sheet = B_sheet_khurana2(orbit_JSO, orbit_SIII_mag)    
    B_ocean = B_induced_infinite(orbit_cphio, B_jupiter + B_sheet, R_C, R_C-surface_layer)
    
    B_total = B_jupiter + B_sheet + B_ocean
    B_total_mag = np.sqrt(B_total[:,0]**2 + B_total[:,1]**2 + B_total[:,2]**2)

    # plot B-field components
    ax_l1.plot(minutes, B_total[:, 0], label="Flyby " + str(flyby_n))
    ax_l2.plot(minutes, B_total[:, 1])
    ax_l3.plot(minutes, B_total[:, 2])
    ax_l4.plot(minutes, B_total_mag)
    
    # plot trajectories
    ax_r1.plot(orbit_cphio[1] / R_C, orbit_cphio[3] / R_C)
    ax_r2.plot(orbit_cphio[1] / R_C, orbit_cphio[2] / R_C)

# add legend
ax_l1.legend(loc="upper right", framealpha=1)

plt.suptitle("Ocean-only Model")
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  

plt.show()
plt.close()

# #-----------plot difference of induction models---------

# # choose flyby
# flyby_n = 6

# orbit_cphio = juice_wrt_callisto["orbit%s" % (flyby_n)]
# orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
# orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
# orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]

# orbit_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')
# orbit_CA = orbit_CA["CA_orbit%s" % (flyby_n)]

# # convert time series to minutes before/after CA
# minutes = [(seconds - orbit_CA[0])/60 for seconds in orbit_cphio[0]]

# #---------compute B-field-----------
# B_jupiter = Bext_Community(orbit_SIII)
# B_sheet = B_sheet_khurana2(orbit_JSO, orbit_SIII_mag)
# B_ocean = B_induced_infinite(orbit_cphio, B_jupiter + B_sheet, R_C, R_C-surface_layer)


# #---------create figure layout-----------

# # Create subplots with 4 rows and 2 columns
# fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(11, 6), sharex=True)

# # clear figure layout
# plt.clf()

# # Specify the width ratios of the columns
# gs = gridspec.GridSpec(4, 2, width_ratios=[3, 1])

# # create subplots for timetraces
# ax_l1 = plt.subplot(gs[0, 0])
# ax_l2 = plt.subplot(gs[1, 0])
# ax_l3 = plt.subplot(gs[2, 0])
# ax_l4 = plt.subplot(gs[3, 0])

# # Create subplots for the trajectories
# ax_r1 = plt.subplot(gs[0:2, 1])
# ax_r2 = plt.subplot(gs[2:4, 1])


# #---------formatting of timetraces-----------

# # add a vertical line at CA
# ax_l1.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
# ax_l2.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
# ax_l3.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)
# ax_l4.axvline(x=0, color='dimgrey', linestyle=":", zorder=0)

# # Apply limits on the time axis
# lim = 45
# time_lim = (-lim, lim)  

# ax_l1.set_xlim(time_lim)
# ax_l2.set_xlim(time_lim)
# ax_l3.set_xlim(time_lim)
# ax_l4.set_xlim(time_lim)

# # set axis labels
# ax_l1.set_ylabel("$B_x$ [nT]")
# ax_l2.set_ylabel("$B_y$ [nT]")
# ax_l3.set_ylabel("$B_z$ [nT]")
# ax_l4.set_ylabel("|B| [nT]")
# ax_l4.set_xlabel("Time after CA [min]")


# #---------formatting of trajectories-----------

# # plot Callisto
# cal1 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)
# cal2 = plt.Circle((0, 0), 1, color='xkcd:dull brown', zorder=1)

# ax_r1.add_patch(cal1)
# ax_r2.add_patch(cal2)

# ax_r1.annotate('Callisto', (1.7, -0.25))
# ax_r2.annotate('Callisto', (1.7, -0.25))

# # Apply limits on the axis
# lim = 6
# orbit_lim = (-lim, lim) 

# ax_r1.set_xlim(orbit_lim)
# ax_r2.set_xlim(orbit_lim)
# ax_r1.set_ylim(orbit_lim)
# ax_r2.set_ylim(orbit_lim)

# # set axis labels
# # ax_r1.set_xlabel('x [$R_C$]')
# ax_r2.set_xlabel('x [$R_C$]')
# ax_r1.set_ylabel('z [$R_C$]')
# ax_r2.set_ylabel('y [$R_C$]')


# #---------general formatting-----------

# # set axis ticks
# ax_l1.tick_params(axis='both', direction='in',top = True, right = True, which='both')
# ax_l2.tick_params(axis='both', direction='in',top = True, right = True, which='both')
# ax_l3.tick_params(axis='both', direction='in',top = True, right = True, which='both')
# ax_l4.tick_params(axis='both', direction='in',top = True, right = True, which='both')

# ax_r1.tick_params(axis='both', direction='in',top = True, right = True, which='both')
# ax_r2.tick_params(axis='both', direction='in',top = True, right = True, which='both')


# ax_l1.xaxis.set_minor_locator(AutoMinorLocator())
# ax_l2.xaxis.set_minor_locator(AutoMinorLocator())
# ax_l3.xaxis.set_minor_locator(AutoMinorLocator())
# ax_l4.xaxis.set_minor_locator(AutoMinorLocator())

# ax_r1.xaxis.set_minor_locator(AutoMinorLocator()) 
# ax_r2.xaxis.set_minor_locator(AutoMinorLocator()) 


# ax_l1.yaxis.set_minor_locator(AutoMinorLocator())
# ax_l2.yaxis.set_minor_locator(AutoMinorLocator())
# ax_l3.yaxis.set_minor_locator(AutoMinorLocator())
# ax_l4.yaxis.set_minor_locator(AutoMinorLocator())

# ax_r1.yaxis.set_minor_locator(AutoMinorLocator()) 
# ax_r2.yaxis.set_minor_locator(AutoMinorLocator()) 


# # plot trajectories
# ax_r1.plot(orbit_cphio[1] / R_C, orbit_cphio[3] / R_C, color="k")
# ax_r2.plot(orbit_cphio[1] / R_C, orbit_cphio[2] / R_C, color="k")

# # define iono conductivities to compare
# iono_conductivities = [0.001, 0.005, 0.01]

# for iono_cond in iono_conductivities:
    
#     conductivities=[0.1, 4, 0.1, iono_cond]
    
#     B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_jupiter + B_sheet, J_omega, conductivities, radii, Styczinski=True)

#     B_diff = B_ocean - B_induced
#     B_diff_mag = np.sqrt(B_diff[:,0]**2 + B_diff[:,1]**2 + B_diff[:,2]**2)

#     # plot B-field components
#     ax_l1.plot(minutes, B_diff[:, 0], label="$\sigma_{iono}$ = %2s mS/m" % (iono_cond*1e3))
#     ax_l2.plot(minutes, B_diff[:, 1])
#     ax_l3.plot(minutes, B_diff[:, 2])
#     ax_l4.plot(minutes, B_diff_mag)
    
# # add legend
# ax_l1.legend(loc="upper right", framealpha=1)

# plt.suptitle("Difference in 4-layer and 3-layer induction model (Flyby %s)" % (flyby_n))
# plt.tight_layout()
# plt.subplots_adjust(hspace=0.25)  

# plt.show()