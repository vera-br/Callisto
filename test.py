# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp

from trajectories.trajectory_analysis import *
from induced_field import *
from jupiter_field import *
from current_sheet import *

# load data
juice_wrt_callisto_cphio = get_spice_data("juice", "callisto", "cphio", "J")
juice_wrt_callisto_cphio_CA = closest_approach_data_J("juice", "callisto", "cphio", "J")
juice_wrt_jupiter_SIII = get_spice_data("juice", "jupiter", "SIII", "J")
#callisto_wrt_jupiter_SIII = get_spice_data("callisto", "jupiter", "SIII", "J")

# specify orbit
orbit_cphio = juice_wrt_callisto_cphio["orbit5"]
orbit_SIII = juice_wrt_jupiter_SIII["orbit5"]
orbit_CA = juice_wrt_callisto_cphio_CA["CA_orbit5"]

# convert time into Timestamp format
OFFSET = datetime(2000,1,1,12) - datetime(1970,1,1) # difference between J2000 and UTC

time = orbit_cphio[0]
time = [Timestamp((datetime.utcfromtimestamp(timestamp) + OFFSET).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]
time_CA = Timestamp((datetime.utcfromtimestamp(orbit_CA[0]) + OFFSET).strftime('%Y-%m-%d %H:%M:%S'))

# radial distance of juice wrt callisto
radial = pd.Series(dict(zip(time, orbit_cphio[4] / R_C)))


#---------magnetic fields-----------

# jovian field
B_external = Bext_full(orbit_SIII)

B_ext_x = pd.Series(dict(zip(time, B_external[:,0])))
B_ext_y = pd.Series(dict(zip(time, B_external[:,1])))
B_ext_z = pd.Series(dict(zip(time, B_external[:,2])))

B_ext_mag = (B_external[:,0]**2 + B_external[:,1]**2 + B_external[:,2]**2)**0.5
B_ext_mag = pd.Series(dict(zip(time, B_ext_mag)))


# induced field
B_induced = B_induced_infinite(orbit_cphio, B_external, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit4_cphio, B_external, 0.03, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)

B_ind_x = pd.Series(dict(zip(time, B_induced[:,0])))
B_ind_y = pd.Series(dict(zip(time, B_induced[:,1])))
B_ind_z = pd.Series(dict(zip(time, B_induced[:,2])))

B_ind_mag = (B_induced[:,0]**2 + B_induced[:,1]**2 + B_induced[:,2]**2)**0.5
B_ind_mag = pd.Series(dict(zip(time, B_ind_mag)))


# current sheet
rho, z, B_sheet = B_disk(orbit_SIII, R0=R0, R1=R1, D=D, I_constant=Icon)

B_sheet_x = pd.Series(dict(zip(time, B_sheet[:,0])))
B_sheet_y = pd.Series(dict(zip(time, B_sheet[:,1])))
B_sheet_z = pd.Series(dict(zip(time, B_sheet[:,2])))

B_sheet_mag = (B_sheet[:,0]**2 + B_sheet[:,1]**2 + B_sheet[:,2]**2)**0.5
B_sheet_mag = pd.Series(dict(zip(time, B_sheet_mag)))

print(B_sheet_mag)





#---------plotting-----------
'''
fig = plt.figure()
ax = fig.gca()
plt.style.use('seaborn-v0_8-whitegrid')

B_ind_mag.plot(ax=ax, label="|B|", color="midnightblue")
B_ind_x.plot(ax=ax, label="Bx", color="royalblue")
B_ind_y.plot(ax=ax, label="By", color="orange")
B_ind_z.plot(ax=ax, label="Bz", color="deeppink")

B_sheet_mag.plot(ax=ax, label="|B|", color="midnightblue")
B_sheet_x.plot(ax=ax, label="Bx", color="royalblue")
B_sheet_y.plot(ax=ax, label="By", color="orange")
B_sheet_z.plot(ax=ax, label="Bz", color="deeppink")

B_ext_mag.plot(ax=ax, label="|B|", color="midnightblue")
B_ext_x.plot(ax=ax, label="Bx", color="royalblue")
B_ext_y.plot(ax=ax, label="By", color="orange")
B_ext_z.plot(ax=ax, label="Bz", color="deeppink")

plt.legend()

radial.plot(ax=ax, label="distance", color="k", secondary_y=True)

# add a vertical line at CA
ax.axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)#, label='closest approach')

# Format the time on the x-axis to include minutes
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))

ax.set_ylabel("B-field [nT]")
ax.right_ax.set_ylabel("Radial distance [R_C]")
ax.set_title("Plasma Sheet field during Flyby 17")
plt.legend(loc="lower right")

plt.show()
'''
