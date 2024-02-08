from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

r_core = 0.09    ; r_surface = 1      ; r_iono = 1.042
sig_core = 1e-6 ; sig_surface = 1e-6 ; sig_iono = 0.5e-3


r_ocean = np.linspace(0.1, 0.99, 100)
sig_ocean = np.linspace( 1e-3, 50e-3, 100)

r_grid, sig_grid = np.meshgrid(r_ocean, sig_ocean)

r_cores = r_core * np.ones_like(r_grid)
r_surfaces = r_surface * np.ones_like(r_grid)
r_ionos = r_iono * np.ones_like(r_grid)
sig_cores = sig_core * np.ones_like(r_grid)
sig_surfaces = sig_surface * np.ones_like(r_grid)
sig_ionos = sig_iono * np.ones_like(r_grid)

conductivities = [sig_cores, sig_grid, sig_surfaces, sig_ionos]
rs = np.array([r_cores, r_grid, r_surfaces, r_ionos]) * R_C

Aeiphi = ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
abs_A = np.abs(Aeiphi)
print(np.shape(abs_A))

phi_A = np.arctan(Aeiphi.real / Aeiphi.imag)
print(np.shape(phi_A))

fig, ax = plt.subplots(1,2)
color1 = ax[0].contourf(r_grid, sig_grid, abs_A)
fig.colorbar(color1)
color2 = ax[1].contourf(r_grid, sig_grid, phi_A)
fig.colorbar(color2)
plt.show()