from induced_field import *
from field_functions import *
import matplotlib.pyplot as plt

# default values
r_core = 0.1 ; r_ocean = 0.9 ; r_surface = 1   ; r_iono = 1.042 
sig_non_conducting = 1e-9 ; sig_ocean = 3 ; sig_iono = 0.5e-3

delta_plot_type = 'conductivity vs conductivity'
conductivity_x = 'Ocean'
conductivity_y = 'Ionosphere'
r_x = 'Ocean'
r_y = 'Ionosphere'

if delta_plot_type == 'conductivity vs radius':
    if conductivity_x == 'Ionosphere':
        sig_x = np.logspace(-5, -1, 100)
    if conductivity_x == 'Ocean':
        sig_x = np.logspace(-3, 1, 100)
    if r_y == 'Ionosphere':
        rs_y = np.linspace(1, 1.2, 100)
    if r_y == 'Ocean':
        rs_y = np.linspace(0.1, 1, 100)
    x_grid, y_grid = np.meshgrid(sig_x, rs_y)

if delta_plot_type == 'conductivity vs conductivity':
    if conductivity_x == 'Ionosphere':
        sig_x = np.logspace(-5, -1, 100)
        sig_y = np.logspace(-3, 1, 100)
    if conductivity_x == 'Ocean':
        sig_x = np.logspace(-3, 1, 100)
        sig_y = np.logspace(-5, -1, 100)
    x_grid, y_grid = np.meshgrid(sig_x, sig_y)

if delta_plot_type == 'radius vs radius':
    if r_x == 'Ionosphere':
        rs_x = np.linspace(1, 1.2, 100)
        rs_y = np.linspace(0.1, 1, 100)
    if r_x == 'Ocean':
        rs_x = np.linspace(0.1, 1, 100)
        rs_y = np.linspace(1, 1.2, 100)
    x_grid, y_grid = np.meshgrid(rs_x, rs_y)



abs_A = np.empty_like(x_grid)
real_A = np.empty_like(x_grid)
phi_A = np.empty_like(x_grid)

abs_A_ocean = np.empty_like(x_grid)
real_A_ocean = np.empty_like(x_grid)
phi_A_ocean = np.empty_like(x_grid)
# conducting ocean and ionosphere

for i in range(np.shape(abs_A)[0]):
    for j in range(np.shape(abs_A)[1]):
        if delta_plot_type == 'conductivity vs radius':
            if conductivity_x == 'Ionosphere':
                conductivities = [sig_non_conducting, sig_ocean, sig_non_conducting, x_grid[i,j]]
                conductivities_ocean = [sig_non_conducting, sig_ocean, sig_non_conducting]
            if conductivity_x == 'Ocean':
                conductivities = [sig_non_conducting, x_grid[i,j], sig_non_conducting, sig_iono]
                conductivities_ocean = [sig_non_conducting, x_grid[i,j], sig_non_conducting]
            if r_y == 'Ionosphere':
                radii = np.array([r_core, r_ocean, r_surface, y_grid[i,j]]) * R_C
                radii_ocean = np.array([r_core, r_ocean, r_surface]) * R_C
            if r_y == 'Ocean':
                radii = np.array([r_core, y_grid[i,j], r_surface, r_iono]) * R_C
                radii_ocean = np.array([r_core, y_grid[i,j], r_surface]) * R_C
        if delta_plot_type == 'conductivity vs conductivity':
            radii = np.array([r_core, r_ocean, r_surface, r_iono]) * R_C
            radii_ocean = np.array([r_core, r_ocean, r_surface]) * R_C
            if conductivity_x == 'Ionosphere':
                conductivities = [sig_non_conducting, y_grid[i,j], sig_non_conducting, x_grid[i,j]]
                conductivities_ocean = [sig_non_conducting, y_grid[i,j], sig_non_conducting]
            if conductivity_x == 'Ocean':
                conductivities = [sig_non_conducting, x_grid[i,j], sig_non_conducting, y_grid[i,j]]
                conductivities_ocean = [sig_non_conducting, x_grid[i,j], sig_non_conducting]
        if delta_plot_type == 'radius vs radius':
            conductivities = [sig_non_conducting, sig_ocean, sig_non_conducting, sig_iono]
            conductivities_ocean = [sig_non_conducting, sig_ocean, sig_non_conducting]
            if r_x == 'Ionosphere':
                radii = np.array([r_core, y_grid[i,j], r_surface, x_grid[i,j]]) * R_C
                radii_ocean = np.array([r_core, y_grid[i,j], r_surface]) * R_C
            if r_x == 'Ocean':
                radii = np.array([r_core, x_grid[i,j], r_surface, y_grid[i,j]]) * R_C
                radii_ocean = np.array([r_core, x_grid[i,j], r_surface]) * R_C

        Aeiphi_both = Aeiphi_Styczinski_many(conductivities, radii, 1, 2*np.pi /(10.1*3600))
        abs_A[i,j] = np.abs(Aeiphi_both)
        real_A[i,j] = Aeiphi_both.real
        phi_A[i,j] = -np.arctan(Aeiphi_both.imag / Aeiphi_both.real) * 180 / np.pi

        Aeiphi_ocean = Aeiphi_Styczinski_many(conductivities_ocean, radii_ocean, 1, 2*np.pi /(10.1*3600))
        abs_A_ocean[i,j] = np.abs(Aeiphi_ocean)
        real_A_ocean[i,j] = Aeiphi_ocean.real
        phi_A_ocean[i,j] = -np.arctan(Aeiphi_ocean.imag / Aeiphi_ocean.real) * 180 / np.pi

# conducting ocean only

delta_abs_A = (abs_A - abs_A_ocean) / abs_A
delta_real_A = (real_A - real_A_ocean) / real_A
delta_phi_A = phi_A - phi_A_ocean

levels_abs_A = np.linspace(0, 1, 11)
levels_phi = np.linspace(-95, 95, 20)
levels_real_A = np.linspace(-1.1, 1.1, 12)

fig, ax = plt.subplots(2,3, layout='constrained')

ax[0,0].contourf(x_grid, y_grid, abs_A, levels_real_A, cmap='seismic')
ax[1,0].contourf(x_grid, y_grid, real_A, levels_real_A, cmap='seismic')

ax[0,1].contourf(x_grid, y_grid, abs_A_ocean, levels_real_A, cmap='seismic')
ax[1,1].contourf(x_grid, y_grid, real_A_ocean, levels_real_A, cmap='seismic')

ax[0,2].contourf(x_grid, y_grid, delta_abs_A, levels_real_A, cmap='seismic')
colorbar = ax[1,2].contourf(x_grid, y_grid, delta_real_A, levels_real_A, cmap='seismic', extend='both')

fig.colorbar(colorbar, ax=ax.ravel().tolist(), ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

fig2, ax2 = plt.subplots(1,3, layout='constrained')

ax2[0].contourf(x_grid, y_grid, phi_A, levels_phi, cmap='hsv')
ax2[1].contourf(x_grid, y_grid, phi_A_ocean, levels_phi, cmap='hsv')
colorbar_phi = ax2[2].contourf(x_grid, y_grid, delta_phi_A, levels_phi, cmap='hsv', extend='both')
fig2.colorbar(colorbar_phi, ax=ax2.ravel().tolist(), ticks=[-90, -60, -30, 0, 30, 60, 90])

ax[0,0].set_title('Ocean and Iono - $|Ae^{i\phi}|$')
ax[1,0].set_title('Ocean and Iono - $Re(Ae^{i\phi})$')
ax[0,1].set_title('Ocean Only - $|Ae^{i\phi}|$')
ax[1,1].set_title('Ocean Only - $Re(Ae^{i\phi})$')
ax[1,2].set_title('$\Delta Re(Ae^{i\phi})$ (%)')
ax[0,2].set_title('$\Delta |Ae^{i\phi}|$ (%)')

ax2[0].set_title('Ocean and Iono - $\phi$')
ax2[1].set_title('Ocean Only - $\phi$')
ax2[2].set_title('$\Delta \phi$')

if delta_plot_type == 'conductivity vs radius':
    for ax_i in ax.ravel():
        ax_i.set_xscale('log')
        ax_i.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
    for ax_i in ax[1]:
        ax_i.set_xlabel('{} Conductivity $\sigma$'.format(conductivity_x))
    for ax_i in ax2:
        ax_i.set_xscale('log')
        ax_i.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
        ax_i.set_xlabel('{} Conductivity $\sigma$'.format(conductivity_x))

if delta_plot_type == 'conductivity vs conductivity':
    for ax_i in ax.ravel():
        ax_i.set_xscale('log')
        ax_i.set_yscale('log')
        ax_i.set_ylabel('{} Conductivity $\sigma$'.format(conductivity_y))
    for ax_i in ax[1]:
        ax_i.set_xlabel('{} Conductivity $\sigma$'.format(conductivity_x))
    for ax_i in ax2:
        ax_i.set_xscale('log')
        ax_i.set_yscale('log')
        ax_i.set_ylabel('{} Conductivity $\sigma$'.format(conductivity_y))
        ax_i.set_xlabel('{} Conductivity $\sigma$'.format(conductivity_x))

if delta_plot_type == 'radius vs radius':
    for ax_i in ax.ravel():
        ax_i.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
    for ax_i in ax[1]:
        ax_i.set_xlabel('{} Outer Radius $[R_C]$'.format(r_x))
    for ax_i in ax2:
        ax_i.set_ylabel('{} Outer Radius $[R_C]$'.format(r_y))
        ax_i.set_xlabel('{} Outer Radius $[R_C]$'.format(r_x))

if delta_plot_type == 'conductivity vs radius':
    fig.suptitle('{} Conductivity vs {} Outer Radius'.format(conductivity_x, r_y))
    fig2.suptitle('{} Conductivity vs {} Outer Radius'.format(conductivity_x, r_y))
if delta_plot_type == 'conductivity vs conductivity':
    fig.suptitle('{} Conductivity vs {} Conductivity'.format(conductivity_x, conductivity_y))
    fig2.suptitle('{} Conductivity vs {} Conductivity'.format(conductivity_x, conductivity_y))
if delta_plot_type == 'radius vs radius':
    fig.suptitle('{} Outer Radius vs {} Outer Radius'.format(r_x, r_y))
    fig2.suptitle('{} Outer Radius vs {} Outer Radius'.format(r_x, r_y))

plt.show()

