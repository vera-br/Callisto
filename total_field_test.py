# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *

# load data

callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII1965')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII_mag')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')

# specify orbit
flyby_n = 2
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
n_interval = 360
phi = np.linspace(0, 2*np.pi, 360)
theta = np.ones_like(phi) * np.pi / 2
r = np.ones_like(phi) * (max(orbit_cal_SIII[4]) + min(orbit_cal_SIII[4])) / 2
orbit_SIII = [theta, theta, theta, theta, r, theta, phi]


#---------magnetic fields-----------

# jovian field
B_external = Bext_Community(orbit_SIII)
Bmag_external = np.sqrt(B_external[:, 0]**2 + B_external[:, 1]**2 + B_external[:, 2]**2)


# current sheet
B_sheet = B_sheet_Community(orbit_SIII)
Bmag_sheet = np.sqrt(B_sheet[:, 0]**2 + B_sheet[:, 1]**2 + B_sheet[:, 2]**2)

B_full_ext = B_external + B_sheet
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)


fig, ax = plt.subplots(2, 2)

ax[0,0].plot(B_external[:,0], label='Jupiter', color='g')
ax[0,0].plot(B_sheet[:,0], label='Sheet', color='m')
ax[0,0].plot(B_full_ext[:,0], label='Full Ext.', color='b')
ax[0,0].set_title('Bx')

ax[0,1].plot(B_external[:,1], label='Jupiter', color='g')
ax[0,1].plot(B_sheet[:,1], label='Sheet', color='m')
ax[0,1].plot(B_full_ext[:,1], label='Full Ext.', color='b')
ax[0,1].set_title('By')


ax[1,0].plot(B_external[:,2], label='Jupiter', color='g')
ax[1,0].plot(B_sheet[:,2], label='Sheet', color='m')
ax[1,0].plot(B_full_ext[:,2], label='Full Ext.', color='b')
ax[1,0].set_title('Bz')

ax[1,1].plot(Bmag_external, label='Jupiter', color='g')
ax[1,1].plot(Bmag_sheet, label='Sheet', color='m')
ax[1,1].plot(Bmag_full_ext, label='Full Ext.', color='b')
ax[1,1].set_title('|B|')
ax[1,1].legend()




plt.show()




