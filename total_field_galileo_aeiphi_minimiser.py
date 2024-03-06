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


# specify orbit
flyby_n = 2
B_polys = []
B_PDSss = []
B_externals = []
orbit_cphios = []

for flyby_n in [1,2]:
    orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
    orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]

    B_PDS = B_PDSs['bfield%s' % (flyby_n)]
    print(np.shape(B_PDS))
    B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

    orbit_cphios.append(orbit_cphio)

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
    B_PDSss.append(np.transpose(B_PDS[1:4]))

# print(np.shape(B_PDSss))
# print(np.shape(B_polys))




def aeiphi_min_func(A, orbit_cphio, B_background, B_PDS):

    measure = 0
    for i in range(2):
        B_induced = B_induced_aeiphi_minimiser(orbit_cphio[i], B_background[i], A)
        B_total_calc = B_background[i] + B_induced
        
        # B_background_norm = (np.array(B_background) - np.array(B_background_mean)) / np.array(B_background_max)
        # B_background_calc_norm = (np.array(B_background_calc) - np.array(B_background_mean[i])) / np.array(B_background_max[i])

        try:
            measure += (root_mean_squared_error(B_PDS[i], B_total_calc))**2
            # measure = 1 - r2_score(B_background, B_background_calc)
        except:
            B_background_calc = np.nan_to_num(B_total_calc, nan=1e30)
            measure += (root_mean_squared_error(B_PDS[i], B_total_calc))**2
            # measure = 1 - r2_score(B_background, B_background_calc)
    return np.sqrt(measure)
# rmse = mean_squared_error(y_test, y_pred, squared = False)
# r2 = r2_score(y_test, y_pred)

from scipy.optimize import minimize
aeiphi_minimised = minimize(aeiphi_min_func, 0.7, args=(orbit_cphios, B_polys, B_PDSss), method='Nelder-Mead')
A = aeiphi_minimised.x

flyby_n = 2
orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]

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

B_induced = B_induced_aeiphi_minimiser(orbit_cphio, B_poly, A)
B_total_calc = B_poly + B_induced

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
ax[0].plot(time, B_total_calc[:,0], color=colour2, linestyle='-', label='Min. Khur.')
ax[1].plot(time, B_total_calc[:,1], color=colour2, linestyle='-')
ax[2].plot(time, B_total_calc[:,2], color=colour2, linestyle='-')

ax[0].set_ylabel('$B_x$ [nT]', fontsize=16)
ax[1].set_ylabel('$B_y$ [nT]', fontsize=16)
ax[2].set_ylabel('$B_z$ [nT]', fontsize=16)

ax[0].set_xlim(min(time), max(time))
ax[1].set_xlim(min(time), max(time))
ax[2].set_xlim(min(time), max(time))

ax[0].legend(framealpha=1, fancybox=True, fontsize=14)

# fig.tight_layout()

plt.show()

