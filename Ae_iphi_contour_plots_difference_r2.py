from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt
from sklearn.metrics import root_mean_squared_error, r2_score
from khurana1997 import *
from jupiter_field import *

# default values

lspace_n = 20

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



# conducting ocean and ionosphere
conductivities = []
radii = []

As = np.linspace(0,1, 10)
phis = np.linspace(0,np.pi/2, 10)

A_grid, phi_grid = np.meshgrid(As, phis)

rsq_err = np.empty_like(A_grid)
r2_grid = np.empty_like(A_grid)
abs_A = np.empty_like(A_grid)
real_A = np.empty_like(A_grid)

for i in range(np.shape(rsq_err)[0]):
    for j in range(np.shape(rsq_err)[1]):
        A_comp = A_grid[i,j]
        aeiphi = A_comp * np.exp(-1j * phi_grid[i,j])
        B_induced = B_induced_finite_conductivity_multilayer_G(orbit_cphio, B_tot_cal, J_omega, conductivities, radii, aeiphi=aeiphi, shifted=True, t_longperiod=t_longperiod)
        rsq = root_mean_squared_error(B_ind, B_induced)
        r2 = r2_score(B_ind, B_induced)

        rsq_err[i,j] = rsq
        r2_grid[i,j] = r2
        abs_A[i,j] = abs(aeiphi)
        real_A[i,j] = aeiphi.real

# conducting ocean only

norm_rsq_err = (rsq_err - rsq_err.min()) / rsq_err.min() * 100
# norm_rsq_err = rsq_err

phi_grid = phi_grid * 180 / np.pi
levels_abs_A = np.linspace(0, 2, 11)
levels_r2 = np.linspace(0, 1, 51)
levels_rsq = np.linspace(0, 100, 21)
levels_real_A = np.linspace(0, 1, 21)

# fig, ax = plt.subplots(1, 3)


# # fig.colorbar(cbar0, ticks=np.linspace(0,1,6))
# cbar1 = ax[0].contourf(A_grid, phi_grid, real_A, levels_real_A, cmap='Reds')
# fig.colorbar(cbar1, ax=ax[0], ticks=np.linspace(0,1,6))
# cs1 = ax[0].contour(A_grid, phi_grid, real_A, levels_real_A, colors='k')

# cbar2 = ax[1].contourf(A_grid, phi_grid, norm_rsq_err, levels_real_A, cmap='Reds_r', extend='max')
# fig.colorbar(cbar2, ticks=np.linspace(0,1,6))
# cs2 = ax[1].contour(A_grid, phi_grid, norm_rsq_err, levels_real_A, colors='k')

# cbar3 = ax[2].contourf(A_grid, phi_grid, r2_grid, levels_real_A, cmap='Reds', extend='min')
# fig.colorbar(cbar3)
# cs3 = ax[2].contour(A_grid, phi_grid, r2_grid, levels_real_A, colors='k')

# ax[0].set_title('$Re(Ae^{-i\phi})$')
# ax[1].set_title('Normalised RMSE')
# ax[2].set_title('$R^2$')


# for ax_i in ax:
#     ax_i.set_xlabel('$|Ae^{-i\phi}|$')
#     ax_i.set_ylabel('$\phi$')

# ax[0].clabel(cs1,manual=True)
# ax[1].clabel(cs2,manual=True)
# ax[2].clabel(cs3,manual=True)

# plt.show()

def fmt(x):
    x = x*100
    s = f"{x:.0f}"
    return f"+{s}%"

fig, ax = plt.subplots(1,1, layout='constrained')
cbar2 = ax.contourf(A_grid, phi_grid, norm_rsq_err, levels_rsq, cmap='Reds_r', extend='max')
cbar = fig.colorbar(cbar2, ticks=np.linspace(0,100,11), label='$\Delta$RMSE (%)')
cbar.ax.set_yticklabels(['+{:.0f}%'.format(x) for x in np.linspace(0,100,11)])

ax.set_xlabel('$|Ae^{-i\phi}|$')
ax.set_ylabel('$\phi$')
plt.show()


