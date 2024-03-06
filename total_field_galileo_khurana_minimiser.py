# load modules and functions
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d 

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from field_functions import *
from jupiter_field import *
from khurana1997 import *

from sklearn.metrics import root_mean_squared_error, r2_score

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

B_polys = []
B_polys_mean = []
B_polys_max = []
B_externals = []
orbit_cal_JSOs = []
orbit_SIII_mags = [] 
orbit_cal_SIIIs = []

for flyby_n in [1,2]:
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

    orbit_cal_JSOs.append(orbit_cal_JSO)
    orbit_SIII_mags.append(orbit_SIII_mag)
    orbit_cal_SIIIs.append(orbit_cal_SIII)

    #----------polynomial fits-------------
    polys = []
    section = int(len(orbit_cphio[0]) / 3)
    for i in range(3):
        t = orbit_cphio[0]
        t = np.append(t[:section], t[2*section:])
        Bi = B_PDS[i+1]
        Bi = np.append(Bi[:section], Bi[2*section:])
        poly = np.polyfit(t, Bi, 3)
        p = np.poly1d(poly)
        polys.append(p)

    B_poly = []
    B_poly_means = []
    B_poly_maxs = []
    for p in polys:
        Bi_poly = p(orbit_cphio[0])
        Bi_poly_mean = np.mean(Bi_poly) * np.ones_like(Bi_poly)
        Bi_poly_max = max(Bi_poly_mean - Bi_poly) * np.ones_like(Bi_poly)
        B_poly_means.append(Bi_poly_mean)
        B_poly_maxs.append(Bi_poly_max)
        B_poly.append(Bi_poly)
    B_poly = np.transpose(B_poly)
    B_poly_mean = np.transpose(B_poly_means)
    B_poly_max = np.transpose(B_poly_maxs)
    B_poly_mag = np.sqrt(B_poly[:, 0]**2 + B_poly[:, 1]**2 + B_poly[:, 2]**2)
    B_polys.append(B_poly)
    B_polys_mean.append(B_poly_mean)
    B_polys_max.append(B_poly_max)

    B_external = Bext_Community(orbit_SIII)
    B_externals.append(B_external)

x0 = -33.5 #* R_J
rho0 = 33.2 #* R_J
v0 = 37.4 #* R_J hr^-1
C1 = 80.3 ; C2 = 690.4 ; C3 = 101.3 ; C4 = -1.7
a1 = 2.49 ; a2 = 1.80 ; a3 = 2.64
r01 = 38.0 #* R_J
rho02 = 2.14 #* R_J
rho03 = 12.5 #* R_J 
D1 = 2.01 #* R_J 
D2 = 13.27 #* R_J
p = 6.26e-3 ; q = 0.35

const_ = [x0, rho0, v0, C1, C2, C3, C4, a1, a2, a3, r01, rho02, rho03, D1, D2, p, q]
bounds_factor = 0.5
bnds = ((1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor), (1-bounds_factor, 1+bounds_factor))
params = np.ones_like(const_)

def khur_min_func(params, orbit_cal_JSO, orbit_SIII_mag, orbit_cal_SIII, B_background, B_background_mean, B_background_max, B_external):
    params = np.multiply(params, const_)
    measure = 0
    for i in range(2):
        B_sheet = B_sheet_khurana(orbit_cal_JSO[i], orbit_SIII_mag[i], orbit_cal_SIII[i], constants=params)
        B_background_calc = B_sheet + B_external[i]
        
        B_background_norm = (np.array(B_background[i]) - np.array(B_background_mean[i])) / np.array(B_background_max[i])
        B_background_calc_norm = (np.array(B_background_calc) - np.array(B_background_mean[i])) / np.array(B_background_max[i])

        try:
            measure += (root_mean_squared_error(B_background_norm, B_background_calc_norm))**2
            # measure = 1 - r2_score(B_background, B_background_calc)
        except:
            B_background_calc = np.nan_to_num(B_background_calc_norm, nan=1e30)
            measure += (root_mean_squared_error(B_background_norm, B_background_calc_norm))**2
            # measure = 1 - r2_score(B_background, B_background_calc)
    return np.sqrt(measure)
# rmse = mean_squared_error(y_test, y_pred, squared = False)
# r2 = r2_score(y_test, y_pred)

from scipy.optimize import minimize
khur_minimised = minimize(khur_min_func, params, args=(orbit_cal_JSOs, orbit_SIII_mags, orbit_cal_SIIIs, B_polys, B_polys_mean, B_polys_max, B_externals), bounds=bnds, method='Nelder-Mead')
constants_min = khur_minimised.x
print(constants_min)
constants_min = constants_min * const_

B_sheet = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag, orbit_cal_SIII, constants=const_)
B_full_ext = B_external + B_sheet
Bmag_full_ext = np.sqrt(B_full_ext[:, 0]**2 + B_full_ext[:, 1]**2 + B_full_ext[:, 2]**2)

rmse_og = root_mean_squared_error(B_poly, B_sheet)

B_sheet_min = B_sheet_khurana(orbit_cal_JSO, orbit_SIII_mag, orbit_cal_SIII, constants=constants_min)
B_full_ext_min = B_external + B_sheet_min
Bmag_full_ext_min = np.sqrt(B_full_ext_min[:, 0]**2 + B_full_ext_min[:, 1]**2 + B_full_ext_min[:, 2]**2)

rmse_min = root_mean_squared_error(B_poly, B_sheet_min)

#---------plot-----------
time = orbit_cphio[0]

fig, ax = plt.subplots(3, 1, sharex=True, constrained_layout=True)

ax[0].plot(time, B_PDS[1], label='Data', color='k', linewidth=0.7, alpha=0.7)
ax[1].plot(time, B_PDS[2], label='PDS', color='k', linewidth=0.7, alpha=0.7)
ax[2].plot(time, B_PDS[3], label='PDS', color='k', linewidth=0.7, alpha=0.7)

colour = 'k'
ax[0].plot(time, B_poly[:,0], color=colour, linestyle='--', label='Poly. Fit')
ax[1].plot(time, B_poly[:,1], color=colour, linestyle='--')
ax[2].plot(time, B_poly[:,2], color=colour, linestyle='--')

colour2 = '#dc267f'
ax[0].plot(time, B_full_ext_min[:,0], color=colour2, linestyle='-', label='Min. Khur.')
ax[1].plot(time, B_full_ext_min[:,1], color=colour2, linestyle='-')
ax[2].plot(time, B_full_ext_min[:,2], color=colour2, linestyle='-')

colour3 = '#785ef0'
ax[0].plot(time, B_full_ext[:,0], color=colour3, linestyle='-', label='OG. Khur.')
ax[1].plot(time, B_full_ext[:,1], color=colour3, linestyle='-')
ax[2].plot(time, B_full_ext[:,2], color=colour3, linestyle='-')

ax[0].set_ylabel('$B_x$ [nT]', fontsize=16)
ax[1].set_ylabel('$B_y$ [nT]', fontsize=16)
ax[2].set_ylabel('$B_z$ [nT]', fontsize=16)

ax[0].set_xlim(min(time), max(time))
ax[1].set_xlim(min(time), max(time))
ax[2].set_xlim(min(time), max(time))

ax[0].legend(framealpha=1, fancybox=True, fontsize=14)

# fig.tight_layout()

plt.show()

