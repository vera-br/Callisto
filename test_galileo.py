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
flyby_n = 4

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()


#callisto

galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
callisto_jupiter_SIII = get_spice_data('callisto', 'jupiter', 'SIII', 'J')
callisto_jupiter_JSO = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')
callisto_jupiter_SIII_mag = get_spice_data('callisto', 'jupiter', 'SIII_mag', 'J')

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)


# create array with time series and constant position in CphiO
time_series = orbit_SIII[0]

constant_pos = np.array([0, 150e3, 0]) #200km above surface on y axis
constant_pos = constant_pos.reshape(1, -1) 
pos_series = np.repeat(constant_pos, len(time_series), axis=0)

orbit_200km = np.c_[time_series, pos_series]

orbit_cphio = orbit_200km.transpose()


#---------magnetic fields-----------

# jovian field
B_external = Bext_Community(orbit_SIII)

# current sheet
B_sheet = B_sheet_khurana2(orbit_JSO, orbit_SIII_mag)


# induced field
r_core = 0.1 * R_C ; r_surface = R_C ; r_iono = 1.05 * R_C; r_ocean = 0.9 * R_C 
sig_core = 1e-9 ; sig_surface = 1e-9 ; sig_iono = 0.01; sig_ocean = 10

radii = [r_core, r_ocean, r_surface, r_iono]

conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]


#B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, 3, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)

B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

# total field
B_total = B_external + B_sheet + B_induced
B_full_ext = B_external + B_sheet


#---------formatting-----------

# calculate B-field magnitudes
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

B_mag_total = np.sqrt(B_total[:,0]**2 + B_total[:,1]**2 + B_total[:,2]**2)
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

# total field
Bx_total = pd.Series(dict(zip(time_model, B_total[:, 0])))
By_total = pd.Series(dict(zip(time_model, B_total[:, 1])))
Bz_total = pd.Series(dict(zip(time_model, B_total[:, 2])))
Bmag_total = pd.Series(dict(zip(time_model, B_mag_total)))

# full external field
Bx_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 0])))
By_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 1])))
Bz_full_ext = pd.Series(dict(zip(time_model, B_full_ext[:, 2])))
Bmag_full_ext = pd.Series(dict(zip(time_model, B_mag_full_ext)))


#---------plotting-----------

fig, ax = plt.subplots(2, 2, figsize=(12,5))

#Bx_Gal.plot(ax=ax[0,0], label="Data", color="skyblue")
Bx_total.plot(ax=ax[0,0], label="External + Induced Field", color="midnightblue")
Bx_full_ext.plot(ax=ax[0,0], label="External Field", color="deeppink")
ax[0,0].set_title('Bx')

#By_Gal.plot(ax=ax[0,1], label="Data", color="skyblue")
By_total.plot(ax=ax[0,1], label="External + Induced Field", color="midnightblue")
By_full_ext.plot(ax=ax[0,1], label="External Field", color="deeppink")
ax[0,1].set_title('By')

#Bz_Gal.plot(ax=ax[1,0], label="Data", color="skyblue")
Bz_total.plot(ax=ax[1,0], label="External + Induced Field", color="midnightblue")
Bz_full_ext.plot(ax=ax[1,0], label="External Field", color="deeppink")
ax[1,0].set_title('Bz')

#Bmag_Gal.plot(ax=ax[1,1], label="Data", color="skyblue")
Bmag_total.plot(ax=ax[1,1], label="External + Induced Field", color="midnightblue")
Bmag_full_ext.plot(ax=ax[1,1], label="External Field", color="deeppink")
ax[1,1].set_title('|B|')

ax[0,0].legend(loc="best", framealpha=1)

# add a vertical line at CA
# ax[0,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[0,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[1,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
# ax[1,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

# set axes labels
ax[0,0].set_ylabel("nT")
ax[0,1].set_ylabel("nT")
ax[1,0].set_ylabel("nT")
ax[1,1].set_ylabel("nT")

# set axes ticks
ax[0,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[0,0].yaxis.set_minor_locator(AutoMinorLocator())

ax[0,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[0,1].yaxis.set_minor_locator(AutoMinorLocator())

ax[1,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1,0].yaxis.set_minor_locator(AutoMinorLocator()) 

ax[1,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1,1].yaxis.set_minor_locator(AutoMinorLocator())

# set grid
ax[0,0].grid(color='xkcd:dark blue',alpha =0.2)
ax[0,1].grid(color='xkcd:dark blue',alpha =0.2)
ax[1,0].grid(color='xkcd:dark blue',alpha =0.2)
ax[1,1].grid(color='xkcd:dark blue',alpha =0.2)

# #title = 'r_iono = ' + str(r_iono / R_C) + ', sig_ocean = ' + str(sig_ocean_2) + ', sig_iono = ' + str(sig_iono)

fig.suptitle("150km above surface")
plt.tight_layout()
plt.show()
