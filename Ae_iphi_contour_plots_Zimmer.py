from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_core = 0.1    ; r_surface = 1    
sig_core = 1e-9

sigma_m = 2 / (mu0 * 2*np.pi /(10.1*3600) * R_C*R_C)

h = np.logspace(-3, 0, 1000)
r_core = np.ones_like(h) - h
sig_ocean = np.logspace(-1, 4, 1000) * sigma_m


# r_core_grid, sig_ocean_grid = np.meshgrid(r_core, sig_ocean)
# r_surfaces = r_surface * np.ones_like(r_core_grid)
# hs = h * np.ones_like(r_core_grid)
# sig_cores = sig_core * np.ones_like(r_core_grid)

sig_ocean_grid, r_core_grid = np.meshgrid(sig_ocean, r_core)
r_surfaces = r_surface * np.ones_like(sig_ocean_grid)
hs = np.ones_like(sig_ocean_grid) - r_core_grid
sig_cores = sig_core * np.ones_like(sig_ocean_grid)

conductivities = [sig_cores, sig_ocean_grid]
rs = np.array([r_core_grid, r_surfaces]) * R_C

Aeiphi = ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
abs_A = np.abs(Aeiphi)
real_A = Aeiphi.real
phi_A = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi

# Aeiphi_abs = np.zeros_like(sig_ocean_grid)
# Aeiphi_phi = np.zeros_like(sig_ocean_grid)
# Aeiphi_real = np.zeros_like(sig_ocean_grid)
# for i in range(np.shape(sig_ocean_grid)[0]):
#     for j in range(np.shape(sig_ocean_grid)[1]):
#         conductivities = [sig_cores[i,j], sig_ocean_grid[i,j]]
#         rs = np.array([r_core_grid[i,j], r_surfaces[i,j]]) * R_C
#         Aeiphi = ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
#         Aeiphi_abs[i,j] = np.abs(Aeiphi)
#         Aeiphi_phi[i,j] = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi
#         Aeiphi_real[i,j] = Aeiphi.real


norm_sig = sig_ocean_grid / sigma_m

fig, ax = plt.subplots(1,3)
color1 = ax[0].contourf(norm_sig, hs, abs_A)
fig.colorbar(color1)
color2 = ax[1].contourf(norm_sig, hs, phi_A)
fig.colorbar(color2)
color3 = ax[2].contourf(norm_sig, hs, real_A)
fig.colorbar(color3)

# fig, ax = plt.subplots(1,3)
# color1 = ax[0].contourf(norm_sig, hs, Aeiphi_abs)
# fig.colorbar(color1)
# color2 = ax[1].contourf(norm_sig, hs, Aeiphi_phi)
# fig.colorbar(color2)
# color3 = ax[2].contourf(norm_sig, hs, Aeiphi_real)
# fig.colorbar(color3)

ax[0].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[0].set_ylabel('Normalised Thickness $r / r_m$')
ax[0].set_title('$|Ae^{i\phi}|$')
ax[0].set_xscale('log')
ax[0].invert_yaxis()
ax[0].set_yscale('log')


ax[1].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[1].set_ylabel('Normalised Thickness $r / r_m$')
ax[1].set_title('$\phi$')
ax[1].set_xscale('log')
ax[1].invert_yaxis()
ax[1].set_yscale('log')

ax[2].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[2].set_ylabel('Normalised Thickness $r / r_m$')
ax[2].set_title('$Re(Ae^{i\phi})$')
ax[2].set_xscale('log')
ax[2].invert_yaxis()
ax[2].set_yscale('log')


#fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
plt.show()

