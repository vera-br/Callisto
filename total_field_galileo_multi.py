# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d
from scipy.optimize import minimize, basinhopping

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()

# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

#----------smoothed B measurements for plotting and r^2-------------------------
Bx_smooth = uniform_filter1d(B_PDS[1], size=500)
By_smooth = uniform_filter1d(B_PDS[2], size=300)
Bz_smooth = uniform_filter1d(B_PDS[3], size=1000)
Bmag_smooth = uniform_filter1d(B_mag, size=300)

Bx_smooth_offset = Bx_smooth - 0.5 * (max(Bx_smooth) + min(Bx_smooth)) * np.ones_like(Bx_smooth)
Bx_smooth_norm = Bx_smooth_offset / max(Bx_smooth_offset)

By_smooth_offset = By_smooth - 0.5 * (max(By_smooth) + min(By_smooth)) * np.ones_like(By_smooth)
By_smooth_norm = By_smooth_offset / max(By_smooth_offset)

Bz_smooth_offset = Bz_smooth - 0.5 * (max(Bz_smooth) + min(Bz_smooth)) * np.ones_like(Bz_smooth)
Bz_smooth_norm = Bz_smooth_offset / max(Bz_smooth_offset)

Bmag_smooth_offset = Bmag_smooth - 0.5 * (max(Bmag_smooth) + min(Bmag_smooth)) * np.ones_like(Bmag_smooth)
Bmag_smooth_norm = Bmag_smooth_offset / max(Bmag_smooth_offset)


#---------magnetic fields-----------

# jovian field
B_external = Bext_Community(orbit_SIII)


# current sheet
B_sheet = B_sheet_Community(orbit_SIII)
B_total = B_external + B_sheet

#-----------multiple induced field parameter run

r_cores = np.array([0.1, 0.2, 0.3]) * R_C ; r_oceans = np.array([0.5, 0.7, 0.9]) * R_C ; r_surface = R_C ; r_ionos = np.array([1.0125, 1.025, 1.0375, 1.05, 1.0625, 1.075]) * R_C
sig_cores = [1e-6]    ; sig_oceans = np.array([0.1, 0.5, 1, 5, 10])      ; sig_surfaces = [1e-6] ; sig_ionos = np.array([0.001, 0.005, 0.01, 0.05, 0.1])


def get_rsq(params, minimise='All', rsq=None, rs=None, sigs=None):
    if minimise == 'All':
        r_core, r_ocean, r_iono, sig_core, sig_ocean, sig_surface, sig_iono = params
    if minimise == 'Radii':
        r_core, r_ocean, r_iono = params
        sig_core, sig_ocean, sig_surface, sig_iono = sigs
    if minimise == 'Conductivities':
        sig_core, sig_ocean, sig_surface, sig_iono = params
        r_core, r_ocean, r_iono = rs

    radii = [r_core, r_ocean, r_surface, r_iono]
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

    # total field
    B_total = B_external + B_sheet + B_induced
    B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

    if rsq == 'rsqx':
        Bx_ind_norm = (B_total[:, 0] - 0.5 * (max(Bx_smooth) + min(Bx_smooth)) * np.ones_like(Bx_smooth)) / max(Bx_smooth_offset)
        Bx_diff = Bx_smooth_norm - Bx_ind_norm
        rsq_ = np.dot(Bx_diff, Bx_diff)
    elif rsq == 'rsqy':
        By_ind_norm = (B_total[:, 1] - 0.5 * (max(By_smooth) + min(By_smooth)) * np.ones_like(By_smooth)) / max(By_smooth_offset)
        By_diff = By_smooth_norm - By_ind_norm
        rsq_ = np.dot(By_diff, By_diff)
    elif rsq == 'rsqz':
        Bz_ind_norm = (B_total[:, 2] - 0.5 * (max(Bz_smooth) + min(Bz_smooth)) * np.ones_like(Bz_smooth)) / max(Bz_smooth_offset)
        Bz_diff = Bz_smooth_norm - Bz_ind_norm
        rsq_ = np.dot(Bz_diff, Bz_diff)
    elif rsq == 'rsqmag':
        Bmag_ind_norm = (B_mag_tot - 0.5 * (max(Bmag_smooth) + min(Bmag_smooth)) * np.ones_like(Bmag_smooth)) / max(Bmag_smooth_offset)
        Bmag_diff = Bmag_smooth_norm - Bmag_ind_norm
        rsq_ = np.dot(Bmag_diff, Bmag_diff)
    elif rsq == 'rsqtot':
        Bx_ind_norm = (B_total[:, 0] - 0.5 * (max(Bx_smooth) + min(Bx_smooth)) * np.ones_like(Bx_smooth)) / max(Bx_smooth_offset)
        By_ind_norm = (B_total[:, 1] - 0.5 * (max(By_smooth) + min(By_smooth)) * np.ones_like(By_smooth)) / max(By_smooth_offset)
        Bz_ind_norm = (B_total[:, 2] - 0.5 * (max(Bz_smooth) + min(Bz_smooth)) * np.ones_like(Bz_smooth)) / max(Bz_smooth_offset)
        Bmag_ind_norm = (B_mag_tot - 0.5 * (max(Bmag_smooth) + min(Bmag_smooth)) * np.ones_like(Bmag_smooth)) / max(Bmag_smooth_offset)
        Bx_diff = Bx_smooth_norm - Bx_ind_norm
        rsqx = np.dot(Bx_diff, Bx_diff)
        By_diff = By_smooth_norm - By_ind_norm
        rsqy = np.dot(By_diff, By_diff)
        Bz_diff = Bz_smooth_norm - Bz_ind_norm
        rsqz = np.dot(Bz_diff, Bz_diff)
        Bmag_diff = Bmag_smooth_norm - Bmag_ind_norm
        rsqmag = np.dot(Bmag_diff, Bmag_diff)
        rsq_ = rsqx + rsqy + rsqz + rsqmag

    return rsq_

# initial guess
r_core_i = 0.15 * R_C ; r_ocean_i = 0.8 * R_C ; r_iono_i = 1.1 * R_C
rs_i = np.array([r_core_i, r_ocean_i, r_iono_i])
sig_core_i = 1e-6 ; sig_ocean_i = 1 ; sig_surface_i = 1e-6 ; sig_iono_i = 0.02
sigs_i = np.array([sig_core_i, sig_ocean_i, sig_surface_i, sig_iono_i])
initial_guess = np.r_[rs_i, sigs_i]
bnds = ((0, R_C), (0, R_C), (R_C, None), (1e-12,None), (1e-12,None), (1e-12,None), (1e-12,None))
bndsx = ((0, R_C), (0, R_C), (R_C, None))
bndssigs = ((1e-12,None), (1e-12,None), (1e-12,None), (1e-12,None))

minimise = 'All'
rsqi = 'rsqx'
tolerance = 1e-9

minimiser_kwargs = {'bounds':bnds}

# note for future: change get_rsq function to take log(sig_i) and to set stepsize

if minimise == 'All':
    resxr = minimize(get_rsq, initial_guess, args=(minimise, rsqi, None, None), bounds=bnds)
    #resxr = basinhopping(get_rsq, initial_guess, minimizer_kwargs={'bounds':bnds, 'args':(minimise, rsqi, None, None)})
    r_core, r_ocean, r_iono, sig_core, sig_ocean, sig_surface, sig_iono = resxr.x
    radii = [r_core, r_ocean, r_surface, r_iono]
    print(np.array(radii)/R_C)
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    print(conductivities)
elif minimise == 'Radii':
    resxr = minimize(get_rsq, rs_i, args=(minimise, rsqi, None, sigs_i), bounds=bndsx)
    r_core, r_ocean, r_iono = resxr.x
    radii = [r_core, r_ocean, r_surface, r_iono]
    print(np.array(radii)/R_C)
    conductivities = sigs_i
elif minimise == 'Conductivities':
    resxsig = minimize(get_rsq, sigs_i, args=(minimise, rsqi, rs_i, None), bounds=bndssigs)
    sig_core, sig_ocean, sig_surface, sig_iono = resxsig.x
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    print(conductivities)
    radii = rs_i

B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

# total field
B_total = B_external + B_sheet + B_induced
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS')
ax[0,1].plot(B_PDS[0], By_smooth, label='PDS')
ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS')
ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS')

ax[0,0].set_title('Bx')
ax[0,1].set_title('By')
ax[1,0].set_title('Bz')
ax[1,1].set_title('|B|')

ax[0,0].plot(B_PDS[0], B_total[:, 0])
ax[0,1].plot(B_PDS[0], B_total[:, 1])
ax[1,0].plot(B_PDS[0], B_total[:, 2])
ax[1,1].plot(B_PDS[0], B_mag_tot, label='Min.')
ax[1,1].legend()

fig.suptitle('r_ocean = ' + str(r_core / R_C) + '-' + str(r_ocean / R_C) + ' R_C, r_iono = 1-' + str(r_iono / R_C) + ' R_C \n  sig_ocean = ' + str(sig_ocean) + ', sig_iono = ' + str(sig_iono))

plt.show()
# plot_all_combs = False

# if plot_all_combs == True:
#     fig, ax = plt.subplots(2, 2)
#     ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS')
#     ax[0,1].plot(B_PDS[0], By_smooth, label='PDS')
#     ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS')
#     ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS')

#     ax[0,0].set_title('Bx')
#     ax[0,1].set_title('By')
#     ax[1,0].set_title('Bz')
#     ax[1,1].set_title('|B|')

# rs = []
# sigs = []
# rsqs = []

# for r_core_i in r_cores:
#     for r_ocean_i in r_oceans:
#         for r_iono_i in r_ionos:
#             for sig_core_i in sig_cores:
#                 for sig_ocean_i in sig_oceans:
#                     for sig_surface_i in sig_surfaces:
#                         for sig_iono_i in sig_ionos:
#                             radii = [r_core_i, r_ocean_i, r_surface, r_iono_i]
#                             conductivities = [sig_core_i, sig_ocean_i, sig_surface_i, sig_iono_i]
#                             B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

#                             # total field
#                             B_total = B_external + B_sheet + B_induced
#                             B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

#                             if plot_all_combs == True:
#                                 ax[0,0].plot(B_PDS[0], B_total[:, 0])
#                                 ax[0,1].plot(B_PDS[0], B_total[:, 1])
#                                 ax[1,0].plot(B_PDS[0], B_total[:, 2])
#                                 ax[1,1].plot(B_PDS[0], B_mag_tot, label=str(r_ocean_i/R_C) + ', ' + str(sig_iono_i))

#                             Bx_ind_norm = (B_total[:, 0] - 0.5 * (max(Bx_smooth) + min(Bx_smooth)) * np.ones_like(Bx_smooth)) / max(Bx_smooth_offset)
#                             By_ind_norm = (B_total[:, 1] - 0.5 * (max(By_smooth) + min(By_smooth)) * np.ones_like(By_smooth)) / max(By_smooth_offset)
#                             Bz_ind_norm = (B_total[:, 2] - 0.5 * (max(Bz_smooth) + min(Bz_smooth)) * np.ones_like(Bz_smooth)) / max(Bz_smooth_offset)
#                             Bmag_ind_norm = (B_mag_tot - 0.5 * (max(Bmag_smooth) + min(Bmag_smooth)) * np.ones_like(Bmag_smooth)) / max(Bmag_smooth_offset)

#                             Bx_diff = Bx_smooth_norm - Bx_ind_norm
#                             rsqx = np.dot(Bx_diff, Bx_diff)
#                             By_diff = By_smooth_norm - By_ind_norm
#                             rsqy = np.dot(By_diff, By_diff)
#                             Bz_diff = Bz_smooth_norm - Bz_ind_norm
#                             rsqz = np.dot(Bz_diff, Bz_diff)
#                             Bmag_diff = Bmag_smooth_norm - Bmag_ind_norm
#                             rsqmag = np.dot(Bmag_diff, Bmag_diff)
#                             rsqtot = rsqx + rsqy + rsqz + rsqmag

#                             rsq_i = [rsqx, rsqy, rsqz, rsqmag, rsqtot]
#                             r_is = [r_core_i, r_ocean_i, r_iono_i]
#                             sig_is = [sig_core_i, sig_ocean_i, sig_surface_i, sig_iono_i]

#                             rsqs.append(rsq_i)
#                             rs.append(r_is)
#                             sigs.append(sig_is)
# if plot_all_combs == True:
#     ax[1,1].legend()

# labels = ['x', 'y', 'z', 'mag', 'tot']

# fig, ax = plt.subplots(2, 2)
# ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS')
# ax[0,1].plot(B_PDS[0], By_smooth, label='PDS')
# ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS')
# ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS')

# ax[0,0].set_title('Bx')
# ax[0,1].set_title('By')
# ax[1,0].set_title('Bz')
# ax[1,1].set_title('|B|')

# for i in range(len(labels)):
#     mindex = np.argmin(np.transpose(rsqs)[i])

#     radii = np.array([rs[mindex][0], rs[mindex][1], r_surface, rs[mindex][2]])
#     conductivities = np.array([sigs[mindex][0], sigs[mindex][1], sigs[mindex][2], sigs[mindex][3]])
#     print(radii / R_C)
#     print(conductivities)

#     B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

#     # total field
#     B_total = B_external + B_sheet + B_induced
#     B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

#     ax[0,0].plot(orbit_cphio[0], B_total[:, 0], label=labels[i])
#     ax[0,1].plot(orbit_cphio[0], B_total[:, 1], label=labels[i])
#     ax[1,0].plot(orbit_cphio[0], B_total[:, 2], label=labels[i])
#     ax[1,1].plot(orbit_cphio[0], B_mag_tot, label=labels[i])

# ax[1,1].legend()
# plt.show()


