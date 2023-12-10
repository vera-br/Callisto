import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp

from trajectories.trajectory_analysis import *
from induced_field import *
from jupiter_field import *


juice_wrt_callisto_cphio = get_spice_data("juice", "callisto", "cphio", "J")
juice_wrt_callisto_cphio_CA = closest_approach_data_J("juice", "callisto", "cphio", "J")
callisto_wrt_jupiter_JSO = get_spice_data("juice", "jupiter", "jupsunorb", "J")

orbit4_cphio = juice_wrt_callisto_cphio["orbit4"]
orbit4_JSO = callisto_wrt_jupiter_JSO["orbit4"]
orbit4_CA = juice_wrt_callisto_cphio_CA["CA_orbit4"]

time = orbit4_cphio[0]
time = [Timestamp(datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]
time_CA = Timestamp(datetime.utcfromtimestamp(orbit4_CA[0]).strftime('%Y-%m-%d %H:%M:%S'))

B_external = Bext_full(orbit4_JSO)
B_induced = B_induced_finite_conductivity(orbit4_cphio, B_external, 0.03, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)

B_ind_x = pd.Series(dict(zip(time, B_induced[:,0])))
B_ind_y = pd.Series(dict(zip(time, B_induced[:,1])))
B_ind_z = pd.Series(dict(zip(time, B_induced[:,2])))

B_ind_mag = (B_induced[:,0]**2 + B_induced[:,1]**2 + B_induced[:,2]**2)**0.5
B_ind_mag = pd.Series(dict(zip(time, B_ind_mag)))

B_ext_x = pd.Series(dict(zip(time, B_external[:,0])))
B_ext_y = pd.Series(dict(zip(time, B_external[:,1])))
B_ext_z = pd.Series(dict(zip(time, B_external[:,2])))

B_ext_mag = (B_external[:,0]**2 + B_external[:,1]**2 + B_external[:,2]**2)**0.5
B_ext_mag = pd.Series(dict(zip(time, B_ext_mag)))

fig = plt.figure()
ax = fig.gca()
plt.style.use('seaborn-v0_8-whitegrid')


B_ind_mag.plot(ax=ax, label="|B|", color="midnightblue")
B_ind_x.plot(ax=ax, label="Bx", color="royalblue")
B_ind_y.plot(ax=ax, label="By", color="orange")
B_ind_z.plot(ax=ax, label="Bz", color="deeppink")
'''
B_ext_mag.plot(ax=ax, label="|B|", color="midnightblue")
B_ext_x.plot(ax=ax, label="Bx", color="royalblue")
B_ext_y.plot(ax=ax, label="By", color="orange")
B_ext_z.plot(ax=ax, label="Bz", color="deeppink")
'''

# Format the time on the x-axis to include minutes
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M'))

# add a vertical line at CA
ax.axvline(x=time_CA, color='dimgrey', linestyle=":", zorder=0)#, label='closest approach')

ax.set_ylabel("B-field [T]")
ax.set_title("Induced B-field during Flyby 4")
plt.legend()
plt.show()

