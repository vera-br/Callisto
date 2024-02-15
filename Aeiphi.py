# calculating the complex amplitude of the induced field

#-------------import modules---------------
import numpy as np
import scipy.special as sps
import scipy.constants as constants
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#-------------import files---------------
from constants import *

#-------------define constants---------------

mu0 = constants.mu_0

R_moon = R_C
surface_layer = 0 #100e3 #m
#ocean_depth = 150e3 #m

sigma_m = 2 / (mu0 * J_omega * R_C**2 )

#-------------define functions---------------

def wavenumber(sigma, omega):
    return (1 - 1j) * np.sqrt(mu0 * sigma * omega / 2)

def zeta(ocean_depth, surface_layer, sigma, omega):
    R_mantle = R_moon - ocean_depth - surface_layer

    z = R_mantle * wavenumber(sigma, omega)

    J_minus_52 = sps.jv(-5 / 2, z)
    J_32 = sps.jv(3 / 2, z)
    J_12 = sps.jv(1 / 2, z)

    return (z * J_minus_52) / (3 * J_32 - z * J_12)

def Aeiphi(ocean_depth, surface_layer, sigma, omega):
    R_ocean = R_moon - surface_layer

    z = R_ocean * wavenumber(sigma, omega)

    J_52 = sps.jv(5 / 2, z)
    J_minus_52 = sps.jv(-5 / 2, z)
    J_12 = sps.jv(1 / 2, z)
    J_minus_12 = sps.jv(-1 / 2, z)

    zet = zeta(ocean_depth, surface_layer, sigma, omega)

    return (R_ocean / R_moon)**3 * (zet*J_52 - J_minus_52) / (zet*J_12 - J_minus_12)


# #-------------testing---------------

# thickness = np.linspace(0, R_moon-surface_layer, 1000)
# conductivities = np.linspace(0, 0.8, 1000)

# CONDUCTIVITY, DEPTH =  np.meshgrid(conductivities, thickness)

# sigma_norm = CONDUCTIVITY / sigma_m
# depth_norm = DEPTH / R_C

# Ae_iphi = Aeiphi(DEPTH, surface_layer, CONDUCTIVITY, J_omega)

# # fig, ax = plt.subplots(1)

# # ax.plot(thickness * 1e-3, -np.real(Ae_iphi))

# fig, ax = plt.subplots()

# ax.set_xlim(0.1, 250)
# ax.set_ylim(0.001, 1)

# plt.xscale('log')
# plt.yscale('log')

# # specify label location
# manual_loc = [(1.8, 0.1), (4.5, 0.07), (9, 0.05), (17, 0.04), (28, 0.033), (43, 0.03), (62, 0.025), (106, 0.02), (160, 0.02)]

# cs = ax.contour(sigma_norm, depth_norm, np.abs(Ae_iphi), levels=np.linspace(0.1,0.9,9))
# ax.clabel(cs, inline=True, fontsize=10, manual=manual_loc)

# ax.invert_yaxis()

# ax.set_xlabel("Normalised Ocean Conductivity")
# ax.set_ylabel("Normalised Ocean Thickness")

# ax.set_title("|Ae^iphi|")

# ax.tick_params(axis='both', direction='in',top = True, right = True, which='both')
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.grid(color='xkcd:dark blue',alpha =0.2)

# plt.show()


# from induced_field import *
# from field_functions import *
# import matplotlib.pyplot as plt

# r_core = 0.1    ; r_surface = 1    
# sig_core = 1e-9

# sigma_m = 2 / (mu0 * 2*np.pi /(10.1*3600) * R_C*R_C)

# h = np.logspace(-3, 0, 1000)
# r_core = np.ones_like(h) - h
# sig_ocean = np.logspace(-1, 4, 1000) * sigma_m


# sig_ocean_grid, r_core_grid = np.meshgrid(sig_ocean, r_core)
# r_surfaces = r_surface * np.ones_like(sig_ocean_grid)
# hs = np.ones_like(sig_ocean_grid) - r_core_grid
# sig_cores = sig_core * np.ones_like(sig_ocean_grid)

# conductivities = [sig_cores, sig_ocean_grid]
# rs = np.array([r_core_grid, r_surfaces]) * R_C

# Aeiphi = #ae_iphi_multilayer(conductivities, rs, 1, 2*np.pi /(10.1*3600))
# abs_A = np.abs(Aeiphi)
# real_A = Aeiphi.real
# phi_A = np.arctan(Aeiphi.real / Aeiphi.imag) * 180 / np.pi


# norm_sig = sig_ocean_grid / sigma_m

# fig, ax = plt.subplots(1,3)
# color1 = ax[0].contourf(norm_sig, hs, abs_A)
# fig.colorbar(color1)
# color2 = ax[1].contourf(norm_sig, hs, phi_A)
# fig.colorbar(color2)
# color3 = ax[2].contourf(norm_sig, hs, real_A)
# fig.colorbar(color3)


# ax[0].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
# ax[0].set_ylabel('Normalised Thickness $r / r_m$')
# ax[0].set_title('$|Ae^{i\phi}|$')
# ax[0].set_xscale('log')
# ax[0].invert_yaxis()
# ax[0].set_yscale('log')


# ax[1].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
# ax[1].set_ylabel('Normalised Thickness $r / r_m$')
# ax[1].set_title('$\phi$')
# ax[1].set_xscale('log')
# ax[1].invert_yaxis()
# ax[1].set_yscale('log')

# ax[2].set_xlabel('Normalised Conductivity $\sigma / \sigma_m$')
# ax[2].set_ylabel('Normalised Thickness $r / r_m$')
# ax[2].set_title('$Re(Ae^{i\phi})$')
# ax[2].set_xscale('log')
# ax[2].invert_yaxis()
# ax[2].set_yscale('log')

# plt.show()
