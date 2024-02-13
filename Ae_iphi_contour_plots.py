from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_core = 0.1    ; r_surface = 1      ; r_iono = 1.042
sig_core = 1e-9 ; sig_surface = 1e-9 ; sig_iono = [0.1e-3, 0.5e-3, 1e-3, 5e-3] ; sig_iono_i = 0.5e-3


r_ocean = np.linspace(0.1, 0.99, 100)
sig_ocean = np.linspace( 1e-3, 50e-3, 100)

r_grid, sig_grid = np.meshgrid(r_ocean, sig_ocean)

r_cores = r_core * np.ones_like(r_grid)
r_surfaces = r_surface * np.ones_like(r_grid)
r_ionos = r_iono * np.ones_like(r_grid)
sig_cores = sig_core * np.ones_like(r_grid)
sig_surfaces = sig_surface * np.ones_like(r_grid)

sig_ionos = sig_iono_i * np.ones_like(r_grid)

conductivities = [sig_cores, sig_grid, sig_surfaces, sig_ionos]
rs = np.array([r_cores, r_grid, r_surfaces, r_ionos]) * R_C

Aeiphi = ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
abs_A = np.abs(Aeiphi)
real_A = Aeiphi.real
phi_A = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi

cons_no_ocean = [sig_cores, sig_ionos]
rs_no_ocean = np.array([r_surfaces, r_ionos]) * R_C

Aeiphi_no_ocean = ae_iphi_multilayer(cons_no_ocean, rs_no_ocean, 1, 2*np.pi /(10.1*3600))
abs_A_no_ocean = np.abs(Aeiphi_no_ocean)
real_A_no_ocean = Aeiphi_no_ocean.real
phi_A_no_ocean = np.arctan(Aeiphi_no_ocean.real / Aeiphi_no_ocean.imag) * 180 / np.pi

fig, ax = plt.subplots(1,3)
color1 = ax[0].contourf(r_grid, sig_grid, abs_A)
fig.colorbar(color1)
color2 = ax[1].contourf(r_grid, sig_grid, phi_A)
fig.colorbar(color2)
color3 = ax[2].contourf(r_grid, sig_grid, real_A)
fig.colorbar(color3)

ax[0].set_xlabel('Ocean radius [$R_C$]')
ax[0].set_ylabel('Ocean conductivity [$Sm^-1$]')
ax[0].set_title('$|Ae^{i\phi}|$')

ax[1].set_xlabel('Ocean radius [$R_C$]')
ax[1].set_ylabel('Ocean conductivity [$Sm^-1$]')
ax[1].set_title('$\phi$')

ax[2].set_xlabel('Ocean radius [$R_C$]')
ax[2].set_ylabel('Ocean conductivity [$Sm^-1$]')
ax[2].set_title('$Re(Ae^{i\phi})$')

fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
plt.show()

# for sig_iono_i in sig_iono:
#     sig_ionos = sig_iono_i * np.ones_like(r_grid)

#     conductivities = [sig_cores, sig_grid, sig_surfaces, sig_ionos]
#     rs = np.array([r_cores, r_grid, r_surfaces, r_ionos]) * R_C

#     Aeiphi = ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
#     abs_A = np.abs(Aeiphi)
#     real_A = Aeiphi.real
#     phi_A = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi

#     cons_no_ocean = [sig_cores, sig_ionos]
#     rs_no_ocean = np.array([r_surfaces, r_ionos]) * R_C

#     Aeiphi_no_ocean = ae_iphi_multilayer(cons_no_ocean, rs_no_ocean, 1, 2*np.pi /(10.1*3600))
#     abs_A_no_ocean = np.abs(Aeiphi_no_ocean)
#     real_A_no_ocean = Aeiphi_no_ocean.real
#     phi_A_no_ocean = np.arctan(Aeiphi_no_ocean.real / Aeiphi_no_ocean.imag) * 180 / np.pi

#     fig, ax = plt.subplots(1,3)
#     color1 = ax[0].contourf(r_grid, sig_grid, abs_A - abs_A_no_ocean)
#     fig.colorbar(color1)
#     color2 = ax[1].contourf(r_grid, sig_grid, phi_A - phi_A_no_ocean)
#     fig.colorbar(color2)
#     color3 = ax[2].contourf(r_grid, sig_grid, real_A - real_A_no_ocean)
#     fig.colorbar(color3)

#     ax[0].set_xlabel('Ocean radius [$R_C$]')
#     ax[0].set_ylabel('Ocean conductivity [$Sm^-1$]')
#     ax[0].set_title('$\Delta|Ae^{i\phi}|$')

#     ax[1].set_xlabel('Ocean radius [$R_C$]')
#     ax[1].set_ylabel('Ocean conductivity [$Sm^-1$]')
#     ax[1].set_title('$\Delta\phi$')

#     ax[2].set_xlabel('Ocean radius [$R_C$]')
#     ax[2].set_ylabel('Ocean conductivity [$Sm^-1$]')
#     ax[2].set_title('$\Delta Re(Ae^{i\phi})$')

#     fig.suptitle('$r_{}$ = {:.1f} $R_C$, $r_{}$ = 1-{:.3f} $R_C$ \n  $\sigma_{}$ = {:2.0e} $Sm^-1$'.format('{ocean}', r_core, '{iono}', r_iono, '{iono}', sig_iono_i))
#     plt.show()