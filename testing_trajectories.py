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
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
galileo_wrt_jupiter_SIII2 = Galileo_trajectories_SIII_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'SIII1965')

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII2["orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)
O_SIII = orbit_SIII.copy()
O_cal_SIII = orbit_cal_SIII.copy()
O_SIII[1:5] = O_SIII[1:5] / R_J
O_cal_SIII[1:5] = O_cal_SIII[1:5] / R_J


fig, ax = plt.subplots(2,3)
ax[0,0].plot(orbit_SIII[0], orbit_SIII[1], 'r')
ax[0,1].plot(orbit_SIII[0], orbit_SIII[2], 'r')
ax[0,2].plot(orbit_SIII[0], orbit_SIII[3], 'r')
ax[1,0].plot(orbit_SIII[0], orbit_SIII[4], 'r')
ax[1,1].plot(orbit_SIII[0], orbit_SIII[5], 'r')
ax[1,2].plot(orbit_SIII[0], orbit_SIII[6], 'r')

ax[0,0].plot(orbit_cal_SIII[0], orbit_cal_SIII[1], 'k')
ax[0,1].plot(orbit_cal_SIII[0], orbit_cal_SIII[2], 'k')
ax[0,2].plot(orbit_cal_SIII[0], orbit_cal_SIII[3], 'k')
ax[1,0].plot(orbit_cal_SIII[0], orbit_cal_SIII[4], 'k')
ax[1,1].plot(orbit_cal_SIII[0], orbit_cal_SIII[5], 'k')
ax[1,2].plot(orbit_cal_SIII[0], orbit_cal_SIII[6], 'k')
plt.show()