from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_core = 0.1    ; r_ocean = 0.9 ; r_surface = 1   ; r_iono = 1.042 
sig_core = 1e-9 ; sig_iono = 0.5e-3

sigma_m = 2 / (mu0 * J_omega * R_C*R_C)

# conducting ocean and ionosphere
# r_ocean = np.linspace(r_core, 0.95, 100)
# sig_ocean = np.logspace(-1, 4, 100) * sigma_m

# sig_grid, r_grid = np.meshgrid(sig_ocean, r_ocean)
# r_cores = r_core * np.ones_like(sig_grid)
# r_surfaces = r_surface * np.ones_like(sig_grid)
# r_ionos = r_iono * np.ones_like(sig_grid)
# sig_cores = sig_core * np.ones_like(sig_grid)
# sig_surfaces = sig_core * np.ones_like(sig_grid)
# sig_ionos = sig_iono * np.ones_like(sig_grid)

# abs_A = np.empty_like(sig_grid)
# real_A = np.empty_like(sig_grid)
# phi_A = np.empty_like(sig_grid)
# for i in range(np.shape(abs_A)[0]):
#     for j in range(np.shape(abs_A)[1]):
#         conductivities = [sig_core, sig_grid[i,j], sig_core, sig_iono]
#         rs = np.array([r_core, r_grid[i,j], r_surface, r_iono]) * R_C
#         Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
#         abs_A[i,j] = np.abs(Aeiphi)
#         real_A[i,j] = Aeiphi.real
#         phi_A[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

# conducting ocean only
r_ocean = np.linspace(r_core, 0.95, 100)
sig_ocean = np.logspace(-1, 4, 100) * sigma_m
sig_grid, r_grid = np.meshgrid(sig_ocean, r_ocean)

abs_A = np.empty_like(sig_grid)
real_A = np.empty_like(sig_grid)
phi_A = np.empty_like(sig_grid)
for i in range(np.shape(abs_A)[0]):
    for j in range(np.shape(abs_A)[1]):
        conductivities = [sig_core, sig_grid[i,j], sig_core]
        rs = np.array([r_core, r_grid[i,j], r_surface]) * R_C
        Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
        abs_A[i,j] = np.abs(Aeiphi)
        real_A[i,j] = Aeiphi.real
        phi_A[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

# ionosphere only
# r_iono = np.linspace(1, 1.1, 100)
# sig_iono = np.logspace(-1, 4, 100) * sigma_m
# sig_grid, r_grid = np.meshgrid(sig_iono, r_iono)

# abs_A = np.empty_like(sig_grid)
# real_A = np.empty_like(sig_grid)
# phi_A = np.empty_like(sig_grid)
# for i in range(np.shape(abs_A)[0]):
#     for j in range(np.shape(abs_A)[1]):
#         conductivities = [sig_core, sig_grid[i,j]]
#         rs = np.array([r_surface, r_grid[i,j]]) * R_C
#         Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
#         abs_A[i,j] = np.abs(Aeiphi)
#         real_A[i,j] = Aeiphi.real
#         phi_A[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

# ionosphere and ocean - conductivities
# r_core = 0.7
# sig_ocean = np.logspace(-2, 2, 100)
# sig_iono = np.logspace(-5, -1, 100)
# sig_grid, sig_grid2 = np.meshgrid(sig_ocean, sig_iono)

# abs_A = np.empty_like(sig_grid)
# real_A = np.empty_like(sig_grid)
# phi_A = np.empty_like(sig_grid)
# for i in range(np.shape(abs_A)[0]):
#     for j in range(np.shape(abs_A)[1]):
#         conductivities = [sig_core, sig_grid[i,j], sig_core, sig_grid2[i,j]]
#         rs = np.array([r_core, r_ocean, r_surface, r_iono]) * R_C
#         Aeiphi = Aeiphi_Styczinski_many(conductivities, rs, 1, 2*np.pi /(10.1*3600))
#         abs_A[i,j] = np.abs(Aeiphi)
#         real_A[i,j] = Aeiphi.real
#         phi_A[i,j] = -np.arctan(Aeiphi.imag / Aeiphi.real) * 180 / np.pi

levels = np.linspace(0,1,11)
levels_real = np.linspace(-1,1, 11)

# x-axis = normalised conductivity, y-axis = outer radius

norm_sig = sig_grid / sigma_m

fig, ax = plt.subplots(1,3)
color1 = ax[0].contourf(norm_sig, r_grid, abs_A, levels)
fig.colorbar(color1)
color2 = ax[1].contourf(norm_sig, r_grid, phi_A)
fig.colorbar(color2)
color3 = ax[2].contourf(norm_sig, r_grid, real_A, levels_real)
fig.colorbar(color3)

ax[0].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[0].set_ylabel('Outer Radius $[R_C]$')
ax[0].set_title('$|Ae^{i\phi}|$')
ax[0].set_xscale('log')

ax[1].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[1].set_ylabel('Outer Radius $[R_C]$')
ax[1].set_title('$\phi$')
ax[1].set_xscale('log')

ax[2].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[2].set_ylabel('Outer Radius $[R_C]$')
ax[2].set_title('$Re(Ae^{i\phi})$')
ax[2].set_xscale('log')

fig, ax = plt.subplots(1,3)
color1 = ax[0].contourf(norm_sig, r_grid, abs_A)
fig.colorbar(color1)
color2 = ax[1].contourf(norm_sig, r_grid, phi_A)
fig.colorbar(color2)
color3 = ax[2].contourf(norm_sig, r_grid, real_A)
fig.colorbar(color3)

ax[0].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[0].set_ylabel('Outer Radius $[R_C]$')
ax[0].set_title('$|Ae^{i\phi}|$')
ax[0].set_xscale('log')

ax[1].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[1].set_ylabel('Outer Radius $[R_C]$')
ax[1].set_title('$\phi$')
ax[1].set_xscale('log')

ax[2].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
ax[2].set_ylabel('Outer Radius $[R_C]$')
ax[2].set_title('$Re(Ae^{i\phi})$')
ax[2].set_xscale('log')

# x-axis - normalised conductivity 1, y-axis - normalised conductivity 2

# sig_grid = sig_grid / sigma_m
# sig_grid2 = sig_grid2 / sigma_m

# fig, ax = plt.subplots(1,3)
# color1 = ax[0].contourf(sig_grid, sig_grid2, abs_A)
# fig.colorbar(color1)
# color2 = ax[1].contourf(sig_grid, sig_grid2, phi_A)
# fig.colorbar(color2)
# color3 = ax[2].contourf(sig_grid, sig_grid2, real_A)
# fig.colorbar(color3)

# ax[0].set_xlabel('Ocean Conductivity $[Sm^{-1}]$')
# ax[0].set_ylabel('Ionosphere Conductivity $[Sm^{-1}]$')
# ax[0].set_title('$|Ae^{i\phi}|$')
# ax[0].set_xscale('log')
# ax[0].set_yscale('log')

# ax[1].set_xlabel('Ocean Conductivity $[Sm^{-1}]$')
# ax[1].set_ylabel('Ionosphere Conductivity $[Sm^{-1}]$')
# ax[1].set_title('$\phi$')
# ax[1].set_xscale('log')
# ax[1].set_yscale('log')

# ax[2].set_xlabel('Ocean Conductivity $[Sm^{-1}]$')
# ax[2].set_ylabel('Ionosphere Conductivity $[Sm^{-1}]$')
# ax[2].set_title('$Re(Ae^{i\phi})$')
# ax[2].set_xscale('log')
# ax[2].set_yscale('log')

#fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
plt.show()

