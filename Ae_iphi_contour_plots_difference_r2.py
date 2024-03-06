from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt
from sklearn.metrics import root_mean_squared_error, r2_score

# default values


lspace_n = 10





flyby_n = 2
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")

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
B_PDS = np.transpose(B_PDS[1:4])
B_ind = B_PDS - B_poly


# conducting ocean and ionosphere
conductivities = []
radii = []

As = np.linspace(0,1, 10)
phis = np.linspace(0,np.pi/2, 10)

A_grid, phi_grid = np.meshgrid(As, phis)

rsq_err = np.empty_like(A_grid)
abs_A = np.empty_like(A_grid)
real_A = np.empty_like(A_grid)

for i in range(np.shape(rsq_err)[0]):
    for j in range(np.shape(rsq_err)[1]):
        A_comp = A_grid[i,j] / np.cos(phi_grid[i,j])
        aeiphi = A_comp * np.exp(-1j * phi_grid[i,j])
        B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_poly, J_omega, conductivities, radii, aeiphi=aeiphi)
        # rsq = root_mean_squared_error(B_ind, B_induced)
        rsq = r2_score(B_ind, B_induced)

        rsq_err[i,j] = rsq
        abs_A[i,j] = abs(aeiphi)
        real_A[i,j] = aeiphi.real

# conducting ocean only

# norm_rsq_err = (rsq_err - rsq_err.min()) / rsq_err.min()
norm_rsq_err = rsq_err

phi_grid = phi_grid * 180 / np.pi
levels_abs_A = np.linspace(0, 2, 11)
levels_rsq = np.linspace(0, 1, 51)
levels_real_A = np.linspace(0, 1, 11)

fig, ax = plt.subplots(1, 3, layout='constrained')


cbar0 = ax[0].contourf(A_grid, phi_grid, abs_A, levels_abs_A, cmap='Reds')
# fig.colorbar(cbar0, ticks=np.linspace(0,1,6))
cbar1 = ax[1].contourf(A_grid, phi_grid, real_A, levels_real_A, cmap='Reds')
fig.colorbar(cbar1, ax=ax[1], ticks=np.linspace(0,1,6))
cbar2 = ax[2].contourf(A_grid, phi_grid, norm_rsq_err, levels_rsq, cmap='seismic', extend='min')
fig.colorbar(cbar2, ticks=np.linspace(0,1,6))


ax[0].set_title('$|Ae^{-i\phi}|$')
ax[1].set_title('$Re(Ae^{-i\phi})$')
ax[2].set_title('$R^2$')


for ax_i in ax:
    ax_i.set_xlabel('A')
    ax_i.set_ylabel('$\phi$')

plt.show()

