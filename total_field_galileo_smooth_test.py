# load modules and functions
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from field_functions import *
from khurana1997 import *
from jupiter_field import *
#from constants import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()


# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]

B_PDS = B_PDSs['bfield%s' % (flyby_n)][:,:7200]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)



J2000 = datetime(2000,1,1,12) # difference between J2000 and UTC
time = orbit_cphio[0][:7200]

# time = [Timestamp((J2000 + timedelta(seconds=timestamp)).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in time]

#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=600)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bx_noise = np.random.normal(size=np.size(Bx_smooth))
By_noise = np.random.normal(size=np.size(By_smooth))
Bz_noise = np.random.normal(size=np.size(Bz_smooth))


#---------plot-----------
# plt.style.use('dark_background')

fig, ax = plt.subplots(3, 1, sharex=True, constrained_layout=True)

ax[0].plot(time, Bx_smooth, label='PDS', color='--r')
ax[1].plot(time, By_smooth, label='PDS', color='--r')
ax[2].plot(time, Bz_smooth, label='PDS', color='--r')

ax[0].plot(time, Bx_smooth + Bx_noise, label='PDS', color='r')
ax[1].plot(time, By_smooth + By_noise, label='PDS', color='r')
ax[2].plot(time, Bz_smooth + Bz_noise, label='PDS', color='r')


ax[0].plot(time, B_PDS[1], label='Data', color='k', linewidth=0.7, alpha=0.7)
ax[1].plot(time, B_PDS[2], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[2].plot(time, B_PDS[3], label='PDS', color='k', linewidth=0.7, alpha=0.7)



ax[0].set_ylabel('$B_x$ [nT]', fontsize=16)
ax[1].set_ylabel('$B_y$ [nT]', fontsize=16)
ax[2].set_ylabel('$B_z$ [nT]', fontsize=16)

# ax[0].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[1].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax[2].xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))

ax[0].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
ax[1].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)
ax[2].tick_params(axis='both', direction='in', top=True, right=True, which='both', labelsize=14)

ax[0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
ax[2].yaxis.set_minor_locator(AutoMinorLocator())

ax[0].set_xlim(min(time), max(time))
ax[1].set_xlim(min(time), max(time))
ax[2].set_xlim(min(time), max(time))

# ax[0].legend(framealpha=1, fancybox=True, fontsize=14)

fig.suptitle('Flyby C9', fontsize=16)

plt.show()
