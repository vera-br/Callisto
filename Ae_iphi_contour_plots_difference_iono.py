from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_core = 0.1    ; r_surface = 1   ; r_iono = 1.042 
sig_core = 1e-9 ; sig_ocean = 3   ; sig_iono = 0.5e-3

# sigma_m = 2 / (mu0 * J_omega * R_C*R_C)

r_ocean = np.linspace(r_core, 0.99, 100)
sig_iono = np.logspace(-5, -1, 1000) #* sigma_m

sig_grid, r_grid = np.meshgrid(sig_iono, r_ocean)
r_cores = r_core * np.ones_like(sig_grid)
r_surfaces = r_surface * np.ones_like(sig_grid)
r_ionos = r_iono * np.ones_like(sig_grid)
sig_cores = sig_core * np.ones_like(sig_grid)
sig_surfaces = sig_core * np.ones_like(sig_grid)
sig_oceans = sig_ocean * np.ones_like(sig_grid)

# conducting ocean and ionosphere

abs_A = np.empty_like(sig_grid)
real_A = np.empty_like(sig_grid)
phi_A = np.empty_like(sig_grid)
for i in range(np.shape(abs_A)[0]):
    for j in range(np.shape(abs_A)[1]):
        conductivities = [sig_core, sig_ocean, sig_core, sig_grid[i,j]]
        rs = np.array([r_core, r_grid[i,j], r_surface, r_iono]) * R_C
        Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
        abs_A[i,j] = np.abs(Aeiphi)
        real_A[i,j] = Aeiphi.real
        phi_A[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

# conducting ocean only

abs_A_ocean = np.empty_like(sig_grid)
real_A_ocean = np.empty_like(sig_grid)
phi_A_ocean = np.empty_like(sig_grid)
for i in range(np.shape(abs_A)[0]):
    for j in range(np.shape(abs_A)[1]):
        conductivities = [sig_core, sig_ocean, sig_core]
        rs = np.array([r_core, r_grid[i,j], r_surface]) * R_C
        Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
        abs_A_ocean[i,j] = np.abs(Aeiphi)
        real_A_ocean[i,j] = Aeiphi.real
        phi_A_ocean[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

delta_abs_A = (abs_A - abs_A_ocean) / abs_A
delta_real_A = (real_A - real_A_ocean) / real_A
delta_phi_A = phi_A - phi_A_ocean

norm_sig = sig_grid # / sigma_m

levels_abs_A = np.linspace(0, 1, 11)
levels_phi = np.linspace(-95, 95, 20)
levels_real_A = np.linspace(-1.1, 1.1, 12)

fig, ax = plt.subplots(2,3)

ax[0,0].contourf(norm_sig, r_grid, abs_A, levels_real_A, cmap='seismic')
ax[1,0].contourf(norm_sig, r_grid, real_A, levels_real_A, cmap='seismic')

ax[0,1].contourf(norm_sig, r_grid, abs_A_ocean, levels_real_A, cmap='seismic')
ax[1,1].contourf(norm_sig, r_grid, real_A_ocean, levels_real_A, cmap='seismic')

ax[0,2].contourf(norm_sig, r_grid, delta_abs_A, levels_real_A, cmap='seismic')
colorbar = ax[1,2].contourf(norm_sig, r_grid, delta_real_A, levels_real_A, cmap='seismic', extend='both')

fig.colorbar(colorbar, ax=ax.ravel().tolist(), ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

ax[0,0].set_ylabel('Ocean Outer Radius $[R_C]$')
ax[0,0].set_title('Ocean and Iono - $|Ae^{i\phi}|$')
ax[0,0].set_xscale('log')

ax[1,0].set_xlabel('Ionosphere Conductivity $\sigma$')
ax[1,0].set_ylabel('Ocean Outer Radius $[R_C]$')
ax[1,0].set_title('Ocean and Iono - $Re(Ae^{i\phi})$')
ax[1,0].set_xscale('log')

ax[0,1].set_ylabel('Ocean Outer Radius $[R_C]$')
ax[0,1].set_title('Ocean Only - $|Ae^{i\phi}|$')
ax[0,1].set_xscale('log')

ax[1,1].set_xlabel('Ionosphere Conductivity $\sigma$')
ax[1,1].set_ylabel('Ocean Outer Radius $[R_C]$')
ax[1,1].set_title('Ocean Only - $Re(Ae^{i\phi})$')
ax[1,1].set_xscale('log')

ax[0,2].set_ylabel('Outer Radius $[R_C]$')
ax[0,2].set_title('$\Delta |Ae^{i\phi}|$ (%)')
ax[0,2].set_xscale('log')

ax[1,2].set_xlabel('Ionosphere Conductivity $\sigma$')
ax[1,2].set_ylabel('Outer Radius $[R_C]$')
ax[1,2].set_title('$\Delta Re(Ae^{i\phi})$ (%)')
ax[1,2].set_xscale('log')

fig2, ax2 = plt.subplots(1,3)

ax2[0].contourf(norm_sig, r_grid, phi_A, levels_phi, cmap='hsv')
ax2[1].contourf(norm_sig, r_grid, phi_A_ocean, levels_phi, cmap='hsv')
colorbar_phi = ax2[2].contourf(norm_sig, r_grid, delta_phi_A, levels_phi, cmap='hsv', extend='both')
fig.colorbar(colorbar_phi, ax=ax2.ravel().tolist(), ticks=[-90, -60, -30, 0, 30, 60, 90])

ax2[0].set_xlabel('Ionosphere Conductivity $\sigma$')
ax2[0].set_ylabel('Ocean Outer Radius $[R_C]$')
ax2[0].set_title('Ocean and Iono - $\phi$')
ax2[0].set_xscale('log')

ax2[1].set_xlabel('Ionosphere Conductivity $\sigma$')
ax2[1].set_ylabel('Ocean Outer Radius $[R_C]$')
ax2[1].set_title('Ocean Only - $\phi$')
ax2[1].set_xscale('log')

ax2[2].set_xlabel('Ionosphere Conductivity $\sigma$')
ax2[2].set_ylabel('Ocean Outer Radius $[R_C]$')
ax2[2].set_title('$\Delta \phi$')
ax2[2].set_xscale('log')


#fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
plt.show()

