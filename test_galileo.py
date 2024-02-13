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

# specify orbit
flyby_n = 2

# load data
# galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
# galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()


# #callisto
# juice_wrt_callisto_cphio = get_spice_data("juice", "callisto", "cphio", "J")
# galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
# callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'G')
# callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'G')
# callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'G')

# orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
# orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
# orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
# orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
# orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
callisto_wrt_jupiter_cphio_CA = get_closest_approach_data("callisto", "jupiter", "SIII", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
galileo_wrt_jupiter_SIII_mag = Galileo_trajectories_SIII_mag_from_CPhiO()
# callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'IAU_JUPITER')
# callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'MAG_VIP4')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = galileo_wrt_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII_CA = callisto_wrt_jupiter_cphio_CA["CA_orbit%s" % (flyby_n)]
# orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
# orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)


# create array with time series and constant position in CphiO
time_series = np.linspace(0, 24, len(orbit_SIII[0]))

constant_pos = np.array([0, 150e3, 0]) #150km above surface on y axis
constant_pos = constant_pos.reshape(1, -1) 
pos_series = np.repeat(constant_pos, len(time_series), axis=0)

orbit_150km = np.c_[time_series, pos_series]

# orbit_cphio = orbit_150km.transpose()


#---------external magnetic field-----------

# jovian field
B_external = Bext_Community(orbit_SIII)

# current sheet
B_sheet = B_sheet_khurana2(orbit_cal_JSO, orbit_SIII_mag)


B_full_ext = B_external + B_sheet


#---------formatting-----------

# calculate B-field magnitudes
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)
B_mag_full_ext = np.sqrt(B_full_ext[:,0]**2 + B_full_ext[:,1]**2 + B_full_ext[:,2]**2)

# convert time from J2000 into Timestamp format
J2000 = datetime(2000,1,1,12)

time_CA = orbit_CA[0]
time_CA = Timestamp((J2000 + timedelta(seconds=orbit_CA[0])).strftime('%Y-%m-%d %H:%M:%S'))

time_model = orbit_cphio[0]
time_model = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time_model]

time_Gal = B_PDS[0]
time_Gal = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time_Gal]

# create series for plotting time evolution
Bx_Gal = pd.Series(dict(zip(time_Gal, B_PDS[1])))
By_Gal = pd.Series(dict(zip(time_Gal, B_PDS[2])))
Bz_Gal = pd.Series(dict(zip(time_Gal, B_PDS[3])))
Bmag_Gal = pd.Series(dict(zip(time_Gal, B_mag)))

# full external field
Bx_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 0])))
By_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 1])))
Bz_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 2])))
Bmag_full_ext = pd.Series(dict(zip(time_model, B_mag_full_ext)))


#---------plotting-----------

fig, ax = plt.subplots(4, 2, figsize=(12,5))

# plot Galileo data and total external field

Bx_Gal.plot(ax=ax[0, 0], label="Data", color="grey")
By_Gal.plot(ax=ax[1, 0], label="Data", color="grey")
Bz_Gal.plot(ax=ax[2, 0], label="Data", color="grey")
Bmag_Gal.plot(ax=ax[3, 0], label="Data", color="grey")

Bx_Gal.plot(ax=ax[0, 1], label="Data", color="grey")
By_Gal.plot(ax=ax[1, 1], label="Data", color="grey")
Bz_Gal.plot(ax=ax[2, 1], label="Data", color="grey")
Bmag_Gal.plot(ax=ax[3, 1], label="Data", color="grey")

Bx_full_ext.plot(ax=ax[0, 0], label="External Field", color="deeppink")
By_full_ext.plot(ax=ax[1, 0], label="External Field", color="deeppink")
Bz_full_ext.plot(ax=ax[2, 0], label="External Field", color="deeppink")
Bmag_full_ext.plot(ax=ax[3, 0], label="External Field", color="deeppink")

Bx_full_ext.plot(ax=ax[0, 1], label="External Field", color="deeppink")
By_full_ext.plot(ax=ax[1, 1], label="External Field", color="deeppink")
Bz_full_ext.plot(ax=ax[2, 1], label="External Field", color="deeppink")
Bmag_full_ext.plot(ax=ax[3, 1], label="External Field", color="deeppink")

# add a vertical line at CA
ax[0, 0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[1, 0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[2, 0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[3, 0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

ax[0, 1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[1, 1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[2, 1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[3, 1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

# set axes labels
ax[0, 0].set_ylabel("Bx (nT)")
ax[1, 0].set_ylabel("By (nT)")
ax[2, 0].set_ylabel("Bz (nT)")
ax[3, 0].set_ylabel("|B| (nT)")

ax[0, 1].set_ylabel("Bx (nT)")
ax[1, 1].set_ylabel("By (nT)")
ax[2, 1].set_ylabel("Bz (nT)")
ax[3, 1].set_ylabel("|B| (nT)")

# set axes ticks
ax[0, 0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1, 0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[2, 0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[3, 0].tick_params(axis='both', direction='in',top = True, right = True, which='both')

ax[0, 1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1, 1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[2, 1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[3, 1].tick_params(axis='both', direction='in',top = True, right = True, which='both')

ax[0, 0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1, 0].yaxis.set_minor_locator(AutoMinorLocator())
ax[2, 0].yaxis.set_minor_locator(AutoMinorLocator())
ax[3, 0].yaxis.set_minor_locator(AutoMinorLocator())

ax[0, 1].yaxis.set_minor_locator(AutoMinorLocator())
ax[1, 1].yaxis.set_minor_locator(AutoMinorLocator())
ax[2, 1].yaxis.set_minor_locator(AutoMinorLocator())
ax[3, 1].yaxis.set_minor_locator(AutoMinorLocator())

# set grid
ax[0, 0].grid(color='xkcd:dark blue',alpha =0.2)
ax[1, 0].grid(color='xkcd:dark blue',alpha =0.2)
ax[2, 0].grid(color='xkcd:dark blue',alpha =0.2)
ax[3, 0].grid(color='xkcd:dark blue',alpha =0.2)

ax[0, 1].grid(color='xkcd:dark blue',alpha =0.2)
ax[1, 1].grid(color='xkcd:dark blue',alpha =0.2)
ax[2, 1].grid(color='xkcd:dark blue',alpha =0.2)
ax[3, 1].grid(color='xkcd:dark blue',alpha =0.2)



#-------------induced field-----------------

# define constants for induction model
surface_layer = 80e3 #m
frequency = 2*np.pi /(10.1*3600)

# vary ocean depth
ocean_depth = np.linspace(0, R_C-surface_layer, 10)
conductivity = 5 # S/m

# set colour map
cmap = plt.colormaps.get_cmap('winter').resampled(len(ocean_depth))
colors = cmap(np.arange(0, cmap.N)) 


for i in range(0, len(ocean_depth)):

    B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, conductivity, frequency, R_C, R_C - surface_layer, R_C -surface_layer - ocean_depth[i])

    B_total = B_external + B_sheet + B_induced
    B_mag_total = np.sqrt(B_total[:,0]**2 + B_total[:,1]**2 + B_total[:,2]**2)

    # combine into series for easy time series formatting
    Bx_total = pd.Series(dict(zip(time_model, B_total[:, 0])))
    By_total = pd.Series(dict(zip(time_model, B_total[:, 1])))
    Bz_total = pd.Series(dict(zip(time_model, B_total[:, 2])))
    Bmag_total = pd.Series(dict(zip(time_model, B_mag_total)))

    Bx_total.plot(ax=ax[0, 0], label=str(ocean_depth[i]), color=colors[i])
    By_total.plot(ax=ax[1, 0], label=str(ocean_depth[i]), color=colors[i])
    Bz_total.plot(ax=ax[2, 0], label=str(ocean_depth[i]), color=colors[i])
    Bmag_total.plot(ax=ax[3, 0], label=str(ocean_depth[i]), color=colors[i])


# vary ocean conductivity
ocean_depth = 150e3 # m
ocean_conductivity = [1,2,3,4,5,6] # S/m

# set colour map
cmap = plt.colormaps.get_cmap('winter').resampled(len(ocean_conductivity))
colors = cmap(np.arange(0, cmap.N)) 


for i in range(0, len(ocean_conductivity)):

    B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, ocean_conductivity[i], frequency, R_C, R_C - surface_layer, R_C -surface_layer - ocean_depth)

    B_total = B_external + B_sheet + B_induced
    B_mag_total = np.sqrt(B_total[:,0]**2 + B_total[:,1]**2 + B_total[:,2]**2)

    # combine into series for easy time series formatting
    Bx_total = pd.Series(dict(zip(time_model, B_total[:, 0])))
    By_total = pd.Series(dict(zip(time_model, B_total[:, 1])))
    Bz_total = pd.Series(dict(zip(time_model, B_total[:, 2])))
    Bmag_total = pd.Series(dict(zip(time_model, B_mag_total)))

    Bx_total.plot(ax=ax[0, 0], label=str(ocean_conductivity[i]), color=colors[i])
    By_total.plot(ax=ax[1, 0], label=str(ocean_conductivity[i]), color=colors[i])
    Bz_total.plot(ax=ax[2, 0], label=str(ocean_conductivity[i]), color=colors[i])
    Bmag_total.plot(ax=ax[3, 0], label=str(ocean_conductivity[i]), color=colors[i])


ax[0, 0].set_title("Varying Ocean depth")
ax[0, 1].set_title("Varying Ocean conductivity")

plt.tight_layout()
plt.show()
