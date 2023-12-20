#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from trajectories.trajectory_analysis import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *

R_C = 2410.3 * 1e3 # m
C_orbit = 1070 * 1e6


# ---------- Coordinate Transformations ---------

# not sure if this works correctly
def spherical_coordinates_transformation(sph_coords, distance):

    r = sph_coords[0]
    theta = sph_coords[1]
    phi = sph_coords[2]

    # Radial Distance
    r_prime = np.sqrt(r**2 + distance**2 - 2 * r * distance * np.cos(theta))

    # Polar Angle
    cos_theta_prime = (r * np.cos(theta) - distance) / r_prime
    sin_theta_prime = (r * np.sin(theta)) / (r_prime * np.sin(phi))

    # Avoid floating-point errors for arccos and arcsin
    cos_theta_prime = np.clip(cos_theta_prime, -1, 1)
    sin_theta_prime = np.clip(sin_theta_prime, -1, 1)

    theta_prime = np.arccos(cos_theta_prime)
    theta_prime = np.where(sin_theta_prime < 0, 2 * np.pi - theta_prime, theta_prime)

    # Azimuthal Angle
    tan_phi_prime = (r * np.sin(theta) * np.sin(phi)) / (r * np.sin(theta) * np.cos(phi) - distance)
    phi_prime = np.arctan(tan_phi_prime)

    return np.column_stack((r_prime, theta_prime, phi_prime))

def spherical_to_cartesian(sph_coords):
    r = sph_coords[:, 0]
    theta = sph_coords[:, 1]
    phi = sph_coords[:, 2]

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return np.column_stack((x, y, z))



# ---------- Magnetic Field ---------

def Bext_Community(x3, y3, z3):
    """
    x,y,z in SIII coords
    """
    jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=False)
    x = x3 / R_J
    y = y3 / R_J
    z = z3 / R_J
    Br, Btheta, Bphi = jm.Internal.Field(x, y, z)
    Bx = Bphi
    By = -Br
    Bz = -Btheta
    return np.array([Bx, By, Bz]).transpose()

def B_sheet_Community(x3, y3, z3):
    """
    x,y,z in SIII coords
    """
    jm.Con2020.Config(equation_type='analytic', CartesianIn=True, CartesianOut=False)
    x = x3 / R_J
    y = y3 / R_J
    z = z3 / R_J
    Br, Btheta, Bphi = jm.Con2020.Field(x, y, z)
    Bx = Bphi
    By = -Br
    Bz = -Btheta
    return np.array([Bx, By, Bz]).transpose()

def B_induced_finite_conductivity_multilayer(position, B_external, omega, conductivities, radii):
    """
    Calculate the induced magnetic field with finite conductivity
    :param position: array with x (m), y(m), z(m) in cphio
    :param Bext_vectors: array of external field vectors Bx, By, Bz in nT
    :param omega: angular frequency of inducing field in
    :param conductivities: array of conductivities of the layers in S
    :param radii: array of radii of the layers in m
    :return: time evolution array of Bx, By, Bz in nT
    """
    A = ae_iphi_multilayer(conductivities, radii, 1, omega).real

    Bind_evolution = []
    for B_ext, pos in zip(B_external, position):

        M = -(2 * pi / mu0) * A * B_ext * (radii[-1]**3)

        rmag = np.linalg.norm(pos)
        rdotM_r = np.dot(pos, M) * pos

        Bind = (mu0 / (4 * pi)) * (3 * rdotM_r - (rmag**2) * M) / (rmag**5)
        Bind_evolution.append(Bind)

    return np.array(Bind_evolution)


# induced field params
r_core = 0.1 * R_C 
r_ocean = 0.9 * R_C 
r_surface = R_C 
r_iono = 1.05 * R_C
sig_core = 1e-9 
sig_ocean = 10 
sig_surface = 1e-9 
sig_iono = 0.01
radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]


# ---------- Setting the Domain ---------

# create time array for one jupiter rotation (synodic)
time = np.arange(0, 10.1*3600, 60)


theta_values = np.linspace(0, np.pi, 100)  # Theta values (latitude)
phi_values = np.linspace(-np.pi, np.pi, 200)  # Phi values (longitude)

# Create a grid of points on the spherical surface
theta, phi = np.meshgrid(theta_values, phi_values)
print(np.shape(theta), np.shape(phi))

B_total_mag = np.zeros((len(phi_values), len(theta_values)))
print(np.shape(B_total_mag))

for j in range(len(theta_values)):
    for i in range(len(phi_values)):

        A = R_C + 50 *1e3

        x = A * np.sin(theta[i,j]) * np.cos(phi[i,j])
        y = A * np.sin(theta[i,j]) * np.sin(phi[i,j])
        z = A * np.cos(theta[i,j])

        position_cphio = np.array([x,y,z])
        sph_cphio = np.array([R_C, theta[i,j], phi[i,j]])

        sph_SIII = spherical_coordinates_transformation(sph_cphio, C_orbit)
        position_SIII = spherical_to_cartesian(sph_SIII)
        position_SIII = position_SIII.transpose()

        B_external = Bext_Community(position_SIII[0], position_SIII[1], position_SIII[2])
        B_sheet = B_sheet_Community(position_SIII[0], position_SIII[1], position_SIII[2])

        B_induced = B_induced_finite_conductivity_multilayer(position_cphio, B_external + B_sheet, omega=2 * np.pi / (10.1 * 3600), conductivities=conductivities, radii=radii)

        B_total = B_external + B_sheet + B_induced
        B_total_mag[i,j] = np.sqrt(B_total[:,0]**2 + B_total[:,1]**2 + B_total[:,2]**2)

print(np.shape(B_total_mag))

#%%
plt.matshow(B_total_mag)
plt.show()
# %%
