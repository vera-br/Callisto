from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_ocean = 0.95 ; r_surface = 1    
sig_core = 1e-9

sigma_m = 2 / (mu0 * 2*np.pi /(10.1*3600) * R_C*R_C)

ocean_depth = np.linspace(0.05, r_ocean, 100) * R_C
r_core = np.ones_like(ocean_depth) * r_ocean * R_C - ocean_depth
#sig_ocean = np.linspace(0.0001, 10, 1000)
sig_ocean = np.logspace(-3, 1, 100)

sig_ocean_grid, ocean_depth_grid = np.meshgrid(sig_ocean, ocean_depth)
surface_layer = (r_surface - r_ocean) * np.ones_like(sig_ocean_grid) * R_C


Aeiphi = Aeiphi(ocean_depth_grid, surface_layer, sig_ocean_grid, J_omega)
abs_A = np.abs(Aeiphi)
real_A = Aeiphi.real
phi_A = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi

levels = np.linspace(0,1,11)
levels_real = np.linspace(-1,1, 11)

fig, ax = plt.subplots(1,3)
color1 = ax[0].contourf(sig_ocean_grid, ocean_depth_grid / R_C, abs_A, levels)
fig.colorbar(color1)
color2 = ax[1].contourf(sig_ocean_grid, ocean_depth_grid / R_C, phi_A)
fig.colorbar(color2)
color3 = ax[2].contourf(sig_ocean_grid, ocean_depth_grid / R_C, real_A, levels_real)
fig.colorbar(color3)

ax[0].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')
ax[1].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')
ax[2].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')

ax[0].set_ylabel('Ocean Thickness $[R_C]$')
ax[1].set_ylabel('Ocean Thickness $[R_C]$')
ax[2].set_ylabel('Ocean Thickness $[R_C]$')

ax[0].set_title('$|Ae^{i\phi}|$')
ax[1].set_title('$\phi$')
ax[2].set_title('$Re(Ae^{i\phi})$')

ax[1].set_xscale('log')
ax[0].set_xscale('log')
ax[2].set_xscale('log')

# fig, ax = plt.subplots(1,3)
# color1 = ax[0].contourf(sig_ocean_grid, ocean_depth_grid / R_C, abs_A)
# fig.colorbar(color1)
# color2 = ax[1].contourf(sig_ocean_grid, ocean_depth_grid / R_C, phi_A)
# fig.colorbar(color2)
# color3 = ax[2].contourf(sig_ocean_grid, ocean_depth_grid / R_C, real_A)
# fig.colorbar(color3)

# ax[0].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')
# ax[1].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')
# ax[2].set_xlabel('Conductivity $\sigma$ $[Sm^{-1}]$')

# ax[0].set_ylabel('Ocean Thickness $[R_C]$')
# ax[1].set_ylabel('Ocean Thickness $[R_C]$')
# ax[2].set_ylabel('Ocean Thickness $[R_C]$')

# ax[0].set_title('$|Ae^{i\phi}|$')
# ax[1].set_title('$\phi$')
# ax[2].set_title('$Re(Ae^{i\phi})$')

#ax[1].set_xscale('log')
#ax[0].set_xscale('log')
#ax[2].set_xscale('log')


#fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
plt.show()

