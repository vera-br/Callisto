# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
from scipy.ndimage import uniform_filter1d
from scipy.optimize import minimize, basinhopping, NonlinearConstraint

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

# # jovian field
# B_external = Bext_Community(orbit_SIII)


# # current sheet
# B_sheet = B_sheet_Community(orbit_SIII)
# B_total = B_external + B_sheet

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

#-----------multiple induced field parameter run
r_surface = R_C ; r_iono = 1.042 * R_C
sig_core = 1e-6 ; sig_surface = 1e-6


def get_rsq(params, minimise='All', rsq=None, rs=None, sigs=None):
    if minimise == 'All':
        r_core, r_ocean = params[0:2] * R_C
        sigs_ = 10**params[2:]
        print(sigs_)
        sig_ocean, sig_iono = sigs_

    if minimise == 'Radii':
        r_core, r_ocean = params * R_C
        sig_ocean, sig_iono = 10**sigs
    if minimise == 'Conductivities':
        sig_ocean, sig_iono = 10**params
        r_core, r_ocean = rs * R_C

    radii = [r_core, r_ocean, r_surface, r_iono]
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, 2*np.pi /(10.1*3600), conductivities, radii)

    # total field
    B_total = B_poly + B_induced
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

#-------------initial guess----------------

r_core_i = 0.3 * R_C ; r_ocean_i = 0.7 * R_C
sig_ocean_i = 1 ; sig_iono_i = 0.1

# __________________________________________

#-------formatting for minimisation---------
rs_i = np.array([r_core_i, r_ocean_i]) / R_C
radii_i = [r_core_i, r_ocean_i, r_surface, r_iono]

sigs_i = np.array([sig_ocean_i, sig_iono_i])
log10_sigs_i = np.log10(sigs_i)

initial_guess = np.r_[rs_i, log10_sigs_i]
bnds = ((0, 1), (0, 1), (-1,2), (-3,1))
bndsx = bnds[0:2]
bndssigs = bnds[2:]

minimise = 'All'
rsqi = 'rsqx'
tolerance = 1e-9

minimiser_kwargs = {'bounds':bnds}

def radii_constraint_func(guess):
    _r_core, _r_ocean = guess[0:2]
    diff = _r_ocean - _r_core
    return diff 
radii_constraint = NonlinearConstraint(radii_constraint_func, 0, np.inf)

if minimise == 'All':
    resxr = minimize(get_rsq, initial_guess, method='trust-constr', constraints=radii_constraint, args=(minimise, rsqi, None, None), bounds=bnds)
    print(resxr.success)
    print(resxr.message)
    
    r_core, r_ocean, log10_sig_ocean, log10_sig_iono = resxr.x
    radii = [r_core * R_C, r_ocean * R_C, r_surface, r_iono]
    sig_ocean, sig_iono = 10**np.array([log10_sig_ocean, log10_sig_iono])
    print(' \n Radii')
    print(np.array(radii) / R_C)
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    print(' \n Conductivities')
    print(conductivities)

elif minimise == 'Radii':
    resxr = minimize(get_rsq, rs_i, args=(minimise, rsqi, None, sigs_i), bounds=bndsx)
    r_core, r_ocean = resxr.x
    radii = [r_core * R_C, r_ocean * R_C, r_surface, r_iono]
    print(np.array(radii)/R_C)
    conductivities = sigs_i

elif minimise == 'Conductivities':
    resxsig = minimize(get_rsq, log10_sigs_i, args=(minimise, rsqi, rs_i, None), bounds=bndssigs)
    sig_core, sig_ocean, sig_surface, sig_iono = resxsig.x
    conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
    print(conductivities)
    radii = rs_i

sigs_init = [sig_core, sig_ocean_i, sig_surface, sig_iono_i]
radii_init = [r_core_i, r_ocean_i, r_surface, r_iono]

B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, 2*np.pi /(10.1*3600), conductivities, radii)
B_induced_init = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, 2*np.pi /(10.1*3600), sigs_init, radii_init)

# total field
B_total = B_poly + B_induced
B_total_i = B_poly + B_induced_init
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)
B_mag_tot_i = np.sqrt(B_total_i[:, 0]**2 + B_total_i[:, 1]**2 + B_total_i[:, 2]**2)

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], Bx_smooth, label='PDS', color='k')
ax[0,1].plot(B_PDS[0], By_smooth, label='PDS', color='k')
ax[1,0].plot(B_PDS[0], Bz_smooth, label='PDS', color='k')
ax[1,1].plot(B_PDS[0], Bmag_smooth, label='PDS Smoothed', color='k')

ax[0,0].plot(B_PDS[0], B_PDS[1], label='PDS', color='k', alpha=0.3)
ax[0,1].plot(B_PDS[0], B_PDS[2], label='PDS', color='k', alpha=0.3)
ax[1,0].plot(B_PDS[0], B_PDS[3], label='PDS', color='k', alpha=0.3)
ax[1,1].plot(B_PDS[0], B_mag, label='PDS', color='k', alpha=0.3)

ax[0,0].set_title('Bx')
ax[0,1].set_title('By')
ax[1,0].set_title('Bz')
ax[1,1].set_title('|B|')

ax[0,0].plot(B_PDS[0], B_poly[:,0], '--r')
ax[0,1].plot(B_PDS[0], B_poly[:,1], '--r')
ax[1,0].plot(B_PDS[0], B_poly[:,2], '--r')
ax[1,1].plot(B_PDS[0], B_poly_mag, '--r')

ax[0,0].plot(B_PDS[0], B_total_i[:, 0], 'b')
ax[0,1].plot(B_PDS[0], B_total_i[:, 1], 'b')
ax[1,0].plot(B_PDS[0], B_total_i[:, 2], 'b')
ax[1,1].plot(B_PDS[0], B_mag_tot_i, 'b', label='Init.')

ax[0,0].plot(B_PDS[0], B_total[:, 0], 'r')
ax[0,1].plot(B_PDS[0], B_total[:, 1], 'r')
ax[1,0].plot(B_PDS[0], B_total[:, 2], 'r')
ax[1,1].plot(B_PDS[0], B_mag_tot, 'r', label='Min.')
ax[1,1].legend()

#fig.suptitle('r_ocean = ' + str(r_core) + '-' + str(r_ocean) + ' R_C, r_iono = 1-' + str(r_iono / R_C) + ' R_C \n  sig_ocean = ' + str(sig_ocean) + ', sig_iono = ' + str(sig_iono))
fig.suptitle('r_ocean = {:.3f} - {:.3f} R_C, r_iono = 1-{:.3f} R_C \n  sig_ocean = {:.3f}, sig_iono = {:.3f}'.format(r_core, r_ocean, r_iono/R_C, sig_ocean, sig_iono))

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


