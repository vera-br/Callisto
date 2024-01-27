# load modules and functions
import matplotlib.pyplot as plt
#from datetime import datetime
#import pandas as pd
#from pandas import Timestamp

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *
from khurana1997 import *


# specify orbit
flyby_n = 2

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
#orbit_SIII_mag = galileo_wrt_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)


#---------magnetic fields-----------

# jovian field
B_external = Bext_Community(orbit_SIII)

# current sheet
#B_sheet = B_sheet_Community(orbit_SIII)
B_sheet = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag)

B_total = B_external + B_sheet

# induced field
r_core = 0.1 * R_C ; r_surface = R_C ; r_iono = 1.05 * R_C
sig_core = 1e-9 ; sig_surface = 1e-9 ; sig_iono = 0.01

r_ocean_1 = 0.8 * R_C 
r_ocean_2 = 0.9 * R_C 
r_ocean_3 = 0.95 * R_C 

radii_1 = [r_core, r_ocean_1, r_surface, r_iono]
radii_2 = [r_core, r_ocean_2, r_surface, r_iono]
radii_3 = [r_core, r_ocean_3, r_surface, r_iono]

sig_ocean_1 = 0.2
sig_ocean_2 = 10
sig_ocean_3 = 500

conductivities_1 = [sig_core, sig_ocean_2, sig_surface, sig_iono]
conductivities_2 = [sig_core, sig_ocean_2, sig_surface, sig_iono]
conductivities_3 = [sig_core, sig_ocean_2, sig_surface, sig_iono]

#B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, 3, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)

B_induced_1 = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities_1, radii_1)
B_induced_2 = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities_2, radii_2)
B_induced_3 = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities_3, radii_3)

# total field
B_total_1 = B_external + B_sheet + B_induced_1
B_total_2 = B_external + B_sheet + B_induced_2
B_total_3 = B_external + B_sheet #+ B_induced_3


#plot_compare_model_with_data(B_total_1, B_PDS, orbit_cphio, orbit_CA, "Galileo Flyby C9")

#---------plots-----------

# calculate B-field magnitudes
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

B_mag_model_1 = np.sqrt(B_total_1[:,0]**2 + B_total_1[:,1]**2 + B_total_1[:,2]**2)
B_mag_model_2 = np.sqrt(B_total_2[:,0]**2 + B_total_2[:,1]**2 + B_total_2[:,2]**2)
B_mag_model_3 = np.sqrt(B_total_3[:,0]**2 + B_total_3[:,1]**2 + B_total_3[:,2]**2)

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

# model 1
Bx_model_1 = pd.Series(dict(zip(time_model, B_total_1[:, 0])))
By_model_1 = pd.Series(dict(zip(time_model, B_total_1[:, 1])))
Bz_model_1 = pd.Series(dict(zip(time_model, B_total_1[:, 2])))
Bmag_model_1 = pd.Series(dict(zip(time_model, B_mag_model_1)))

# model 2
Bx_model_2 = pd.Series(dict(zip(time_model, B_total_2[:, 0])))
By_model_2 = pd.Series(dict(zip(time_model, B_total_2[:, 1])))
Bz_model_2 = pd.Series(dict(zip(time_model, B_total_2[:, 2])))
Bmag_model_2 = pd.Series(dict(zip(time_model, B_mag_model_2)))

# model 3
Bx_model_3 = pd.Series(dict(zip(time_model, B_total_3[:, 0])))
By_model_3 = pd.Series(dict(zip(time_model, B_total_3[:, 1])))
Bz_model_3 = pd.Series(dict(zip(time_model, B_total_3[:, 2])))
Bmag_model_3 = pd.Series(dict(zip(time_model, B_mag_model_3)))

# plot
fig, ax = plt.subplots(2, 2, figsize=(12,5))

Bx_Gal.plot(ax=ax[0,0], label="Data", color="skyblue")
Bx_model_1.plot(ax=ax[0,0], label="External + Induced Field", color="midnightblue")
Bx_model_3.plot(ax=ax[0,0], label="External Field", color="deeppink")
ax[0,0].set_title('Bx')

By_Gal.plot(ax=ax[0,1], label="Data", color="skyblue")
By_model_1.plot(ax=ax[0,1], label="External + Induced Field", color="midnightblue")
#By_model_2.plot(ax=ax[0,1], label="surface = " + str(0.1 * R_C / 1e3) + "km", color="gold", alpha=.3)
By_model_3.plot(ax=ax[0,1], label="External Field", color="deeppink")
ax[0,1].set_title('By')

Bz_Gal.plot(ax=ax[1,0], label="Data", color="skyblue")
Bz_model_1.plot(ax=ax[1,0], label="External + Induced Field", color="midnightblue")
#Bz_model_2.plot(ax=ax[1,0], label="surface = " + str(0.1 * R_C / 1e3) + "km", color="gold", alpha=.3)
Bz_model_3.plot(ax=ax[1,0], label="External Field", color="deeppink")
ax[1,0].set_title('Bz')

Bmag_Gal.plot(ax=ax[1,1], label="Data", color="skyblue")
Bmag_model_1.plot(ax=ax[1,1], label="External + Induced Field", color="midnightblue")
#Bmag_model_2.plot(ax=ax[1,1], label="surface = " + str(0.1 * R_C / 1e3) + "km", color="gold", alpha=.3)
Bmag_model_3.plot(ax=ax[1,1], label="External Field", color="deeppink")
ax[1,1].set_title('|B|')

ax[0,0].legend(loc="best", framealpha=1)

# add a vertical line at CA
ax[0,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[0,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[1,0].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)
ax[1,1].axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)

# set axes labels
ax[0,0].set_ylabel("nT")
ax[0,1].set_ylabel("nT")
ax[1,0].set_ylabel("nT")
ax[1,1].set_ylabel("nT")

# set axes ticks
ax[0,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[0,0].yaxis.set_minor_locator(AutoMinorLocator())
#ax[0,0].xaxis.set_minor_locator(AutoMinorLocator())

ax[0,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[0,1].yaxis.set_minor_locator(AutoMinorLocator())
#ax[0,1].xaxis.set_minor_locator(AutoMinorLocator())

ax[1,0].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1,0].yaxis.set_minor_locator(AutoMinorLocator()) 
#ax[1,0].xaxis.set_minor_locator(AutoMinorLocator())

ax[1,1].tick_params(axis='both', direction='in',top = True, right = True, which='both')
ax[1,1].yaxis.set_minor_locator(AutoMinorLocator())
#ax[1,1].xaxis.set_minor_locator(AutoMinorLocator())

# set grid
ax[0,0].grid(color='xkcd:dark blue',alpha =0.2)
ax[0,1].grid(color='xkcd:dark blue',alpha =0.2)
ax[1,0].grid(color='xkcd:dark blue',alpha =0.2)
ax[1,1].grid(color='xkcd:dark blue',alpha =0.2)

#title = 'r_iono = ' + str(r_iono / R_C) + ', sig_ocean = ' + str(sig_ocean_2) + ', sig_iono = ' + str(sig_iono)

#fig.suptitle(title)
plt.tight_layout()
plt.show()
