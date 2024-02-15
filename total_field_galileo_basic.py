# load modules and functions
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from field_functions import *
#from constants import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
callisto_wrt_jupiter_cphio_CA = get_closest_approach_data("callisto", "jupiter", "SIII", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()
galileo_wrt_jupiter_SIII_mag = Galileo_trajectories_SIII_mag_from_CPhiO()
callisto_jupiter_SIII = find_nearest_trajectories_G('callisto', 'jupiter', 'IAU_JUPITER')
callisto_jupiter_SIII_mag = find_nearest_trajectories_G('callisto', 'jupiter', 'MAG_VIP4')
callisto_jupiter_JSO = find_nearest_trajectories_G('callisto', 'jupiter', 'jupsunorb')

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = galileo_wrt_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII_CA = callisto_wrt_jupiter_cphio_CA["CA_orbit%s" % (flyby_n)]
orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

#----------smoothed B measurements for plotting-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

#----------polynomial fits-------------
polys = []
section = int(len(B_PDS[0]) / 3)
for i in range(3):
    t = B_PDS[0]
    t = np.append(t[:section], t[2*section:])
    Bi = B_PDS[i+1]
    Bi = np.append(Bi[:section], Bi[2*section:])
    poly = np.polyfit(t, Bi, 3)
    p = np.poly1d(poly)
    polys.append(p)

B_poly = []
for p in polys:
    Bi_poly = p(B_PDS[0])
    B_poly.append(Bi_poly)
B_poly = np.transpose(B_poly)
B_poly_mag = np.sqrt(B_poly[:, 0]**2 + B_poly[:, 1]**2 + B_poly[:, 2]**2)


# induced field parameters

r_core = 0.7 * R_C ;   r_ocean = 0.9 * R_C ;   r_surface = R_C
sig_ocean = 1e0

ocean_depth = r_ocean - r_core
surface_layer = r_surface - r_ocean


# Induced Field calculation
# B_induced_fincon = B_induced_finite_conductivity2(orbit_cphio, B_poly, sig_ocean, J_omega, 0.2*R_C, 0.1*R_C)
B_induced_fincon = B_induced_infinite(orbit_cphio, B_poly, 0.7*R_C, R_C)
Bmag_induced_fincon = np.sqrt(B_induced_fincon[:, 0]**2 + B_induced_fincon[:, 1]**2 + B_induced_fincon[:, 2]**2)
B_total_fincon = B_poly - B_induced_fincon
B_mag_tot_fincon = np.sqrt(B_total_fincon[:, 0]**2 + B_total_fincon[:, 1]**2 + B_total_fincon[:, 2]**2)

#---------plot-----------

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS Smoothed', color='k')

ax[0,0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax[0,1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
ax[1,0].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
ax[1,1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)

# Polynomial
ax[0,0].plot(B_PDS[0], B_poly[:,0], '--r')
ax[0,1].plot(B_PDS[0], B_poly[:,1], '--r')
ax[1,0].plot(B_PDS[0], B_poly[:,2], '--r')
ax[1,1].plot(B_PDS[0], B_poly_mag, '--r', label='External')

ax[0,0].plot(orbit_cphio[0], B_total_fincon[:, 0], color='r')
ax[0,1].plot(orbit_cphio[0], B_total_fincon[:, 1], color='r')
ax[1,0].plot(orbit_cphio[0], B_total_fincon[:, 2], color='r')
ax[1,1].plot(orbit_cphio[0], B_mag_tot_fincon, label='Total', color='r')

ax[0,0].set_title('$B_x$')
ax[0,1].set_title('$B_y$')
ax[1,0].set_title('$B_z$')
ax[1,1].set_title('$|B|$')

ax[0,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[0,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[1,0].set_xlim(min(B_PDS[0]), max(B_PDS[0]))
ax[1,1].set_xlim(min(B_PDS[0]), max(B_PDS[0]))

ax[1,1].legend(framealpha=1, fancybox=True)


fig.suptitle('$r_{} = {:.1f}-{:.1f} R_C, \sigma_{} = {:2.0e} Sm^{}$'.format('{ocean}', r_core / R_C, r_ocean / R_C, '{ocean}', sig_ocean, '{-1}'))
plt.show()

