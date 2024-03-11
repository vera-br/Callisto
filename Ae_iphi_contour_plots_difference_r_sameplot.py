from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt
from sklearn.metrics import root_mean_squared_error, r2_score
from khurana1997 import *
from jupiter_field import *
import matplotlib.ticker as mtick 

# default values
r_core = 0.8 ; r_ocean = 0.9 ; r_surface = 1   ; r_iono = 1.042 
sig_non_conducting = 1e-9 ; sig_ocean = 3 ; sig_iono = 0.5e-5

delta_plot_type = 'conductivity vs radius'
conductivity_x = 'Ocean'
r_x = 'Ocean'

conductivity_y = 'Ionosphere'
r_y = 'Ocean'

lspace_n = 20

sig_oceans = np.logspace(-3, 2, lspace_n)
sig_ionos = np.logspace(-5, -1, lspace_n)
rs_oceans = np.linspace(r_core, 1, lspace_n)
rs_ionos = np.linspace(1, 1.2, lspace_n)

if delta_plot_type == 'conductivity vs radius':
    if conductivity_x == 'Ionosphere':
        sig_x = sig_ionos
    if conductivity_x == 'Ocean':
        sig_x = sig_oceans
    if r_y == 'Ionosphere':
        rs_y = rs_ionos
    if r_y == 'Ocean':
        rs_y = rs_oceans
    x_grid, y_grid = np.meshgrid(sig_x, rs_y)

if delta_plot_type == 'conductivity vs conductivity':
    if conductivity_x == 'Ionosphere':
        sig_x = sig_ionos
        sig_y = sig_oceans
    if conductivity_x == 'Ocean':
        sig_x = sig_oceans
        sig_y = sig_ionos
    x_grid, y_grid = np.meshgrid(sig_x, sig_y)

if delta_plot_type == 'radius vs radius':
    if r_x == 'Ionosphere':
        rs_x = rs_ionos
        rs_y = rs_oceans
    if r_x == 'Ocean':
        rs_x = rs_oceans
        rs_y = rs_ionos
    x_grid, y_grid = np.meshgrid(rs_x, rs_y)



flyby_n = 2
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")

callisto_jupiter_SIII_longperiod = get_spice_data('callisto', 'jupiter', 'SIII_fullcycle', 'G')
callisto_jupiter_SIII_mag_longperiod = get_spice_data('callisto', 'jupiter', 'SIII_mag_fullcycle', 'G')
callisto_jupiter_JSO_longperiod = get_spice_data('callisto', 'jupiter', 'jupsunorb_fullcycle', 'G')

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_cal_SIII_LP = callisto_jupiter_SIII_longperiod["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag_LP = callisto_jupiter_SIII_mag_longperiod["orbit%s" % (flyby_n)]
orbit_cal_JSO_LP = callisto_jupiter_JSO_longperiod["orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]

t_longperiod = orbit_cal_SIII_LP[0]

B_jup = Bext_Community(orbit_cal_SIII_LP)
B_sheet = B_sheet_khurana(orbit_cal_JSO_LP, orbit_cal_SIII_mag_LP, orbit_cal_SIII_LP)
B_tot_cal = B_jup + B_sheet

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
B_PDS = np.transpose(B_PDS[1:4])
B_ind = B_PDS - B_poly

rsq_err = np.empty_like(x_grid)
r2_grid = np.empty_like(x_grid)
abs_A = np.empty_like(x_grid)
real_A = np.empty_like(x_grid)
phi_grid = np.empty_like(x_grid)
# conducting ocean and ionosphere

for i in range(np.shape(rsq_err)[0]):
    for j in range(np.shape(rsq_err)[1]):
        if delta_plot_type == 'conductivity vs radius':
            if conductivity_x == 'Ionosphere':
                conductivities = [sig_non_conducting, sig_ocean, sig_non_conducting, x_grid[i,j]]

            if conductivity_x == 'Ocean':
                conductivities = [sig_non_conducting, x_grid[i,j], sig_non_conducting, sig_iono]
             
            if r_y == 'Ionosphere':
                radii = np.array([r_core, r_ocean, r_surface, y_grid[i,j]]) * R_C
               
            if r_y == 'Ocean':
                radii = np.array([r_core, y_grid[i,j], r_surface, r_iono]) * R_C
                
        if delta_plot_type == 'conductivity vs conductivity':
            radii = np.array([r_core, r_ocean, r_surface, r_iono]) * R_C
            
            if conductivity_x == 'Ionosphere':
                conductivities = [sig_non_conducting, y_grid[i,j], sig_non_conducting, x_grid[i,j]]
               
            if conductivity_x == 'Ocean':
                conductivities = [sig_non_conducting, x_grid[i,j], sig_non_conducting, y_grid[i,j]]
             
        if delta_plot_type == 'radius vs radius':
            conductivities = [sig_non_conducting, sig_ocean, sig_non_conducting, sig_iono]
         
            if r_x == 'Ionosphere':
                radii = np.array([r_core, y_grid[i,j], r_surface, x_grid[i,j]]) * R_C
          
            if r_x == 'Ocean':
                radii = np.array([r_core, x_grid[i,j], r_surface, y_grid[i,j]]) * R_C
               
        aeiphi = Aeiphi_Styczinski_many(conductivities, radii, 1, J_omega)
        B_induced = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_tot_cal, J_omega, conductivities, radii, aeiphi=aeiphi, shifted=True, t_longperiod=t_longperiod)
        rsq = root_mean_squared_error(B_ind, B_induced)
        r2 = r2_score(B_ind, B_induced)

        rsq_err[i,j] = rsq
        r2_grid[i,j] = r2
        abs_A[i,j] = abs(aeiphi)
        real_A[i,j] = aeiphi.real
        phi_grid[i,j] = -np.angle(aeiphi)

# conducting ocean only

norm_rsq_err = (rsq_err - rsq_err.min()) / rsq_err.min()
# norm_rsq_err = rsq_err

phi_grid = phi_grid * 180 / np.pi
levels_abs_A = np.linspace(0, 2, 11)

levels_r2 = np.linspace(0, 1, 51)
levels_r2_contour = np.linspace(0, 1, 11)

levels_rsq = np.linspace(0, 1, 51)
levels_rsq_contour = np.linspace(0, 1, 11)

levels_real_A = np.linspace(0, 1, 11)
levels_real_A_contour = np.linspace(0, 1, 11)

levels_phi = np.linspace(0, 90, 19)
levels_phi_contour = np.linspace(-95, 95, 10)

fig, ax = plt.subplots(1, 2, layout='constrained')

cbar0 = ax[0].contourf(x_grid, y_grid, abs_A, levels_real_A, cmap='Reds')
cs0 = ax[0].contour(x_grid, y_grid, norm_rsq_err, levels_rsq_contour, colors='k')
fig.colorbar(cbar0, ax=ax[0], orientation='horizontal', ticks=np.linspace(0,1,6), label='$|A|$')

cbar2 = ax[1].contourf(x_grid, y_grid, phi_grid, levels_phi, cmap='Reds')
cs1 = ax[1].contour(x_grid, y_grid, norm_rsq_err, levels_rsq_contour, colors='k')
fig.colorbar(cbar2, ax=ax[1], orientation='horizontal', ticks=np.linspace(0, 90, 10), label='$\phi$')

# ax[0].set_title('$|Ae^{-i\phi}|$')
# ax[1].set_title('$Re(Ae^{-i\phi})$')
# ax[2].set_title('$\phi$')
# ax[3].set_title('Normalised RMSE')
# ax[4].set_title('$R^2$')

if delta_plot_type == 'conductivity vs radius':
    for ax in ax.ravel():
        ax.set_xscale('log')
        ax.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
        ax.set_xlabel('{} Conductivity $\sigma$ $[Sm^{}]$'.format(conductivity_x, -1))


if delta_plot_type == 'conductivity vs conductivity':
    for ax in ax.ravel():
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel('{} Conductivity $\sigma$ $[Sm^{}]$'.format(conductivity_y))
        ax.set_xlabel('{} Conductivity $\sigma$ $[Sm^{}]$'.format(conductivity_x))


if delta_plot_type == 'radius vs radius':
    for ax in ax.ravel():
        ax.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
        ax.set_xlabel('{} Outer Radius $[R_C]$'.format(r_x))


if delta_plot_type == 'conductivity vs radius':
    fig.suptitle('{} Conductivity vs {} Outer Radius'.format(conductivity_x, r_y))
if delta_plot_type == 'conductivity vs conductivity':
    fig.suptitle('{} Conductivity vs {} Conductivity'.format(conductivity_x, conductivity_y))
if delta_plot_type == 'radius vs radius':
    fig.suptitle('{} Outer Radius vs {} Outer Radius'.format(r_x, r_y))

def fmt(x):
    x = x*100
    s = f"{x:.0f}"
    return f"+{s}%"

plt.clabel(cs0, cs0.levels, manual=True, fmt=fmt)
plt.clabel(cs1, cs1.levels, manual=True, fmt=fmt)

plt.show()

