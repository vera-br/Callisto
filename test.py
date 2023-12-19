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

# ---------- Setting the Domain ---------

# define points on sphere surface
def points_on_sphere(num_points, radius=R_C):
    # Generate random spherical coordinates
    theta = np.random.uniform(0, 2 * np.pi, num_points)
    phi = np.arccos(2 * np.random.uniform(0, 1, num_points) - 1)

    # Convert spherical coordinates to Cartesian coordinates
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    # Return the points as a NumPy array
    return x, y, z

# create time array for one jupiter rotation (synodic)
time = np.arange(0, 10.1*3600, 60)

x, y, z = points_on_sphere(1000, R_C)

domain_cphio = np.empty((1,4))
for i in range(len(time)):
    t = np.full((len(x)), fill_value=time[i])
    domain_cphio = np.concatenate((domain_cphio, np.column_stack((t,x,y,z))))
domain_cphio = np.delete(domain_cphio, (0), axis=0) #delete first row (created with np.empty)

time_array = domain_cphio[:,0]

# ---------- Coordinate Transformations ---------

# convert carttesian to sph (cphio)
sph_cphio = cartesian_to_spherical(domain_cphio[:,1:4])
sph_domain_calisto = np.concatenate((time_array.reshape(-1, 1), sph_cphio), axis=1)

# not sure if this works correctly
def spherical_coordinates_transformation(sph_coords, distance):
    r = sph_coords[:, 0]
    theta = sph_coords[:, 1]
    phi = sph_coords[:, 2]

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

# transform from sph_cphio to sph_SIII
sph_SIII = spherical_coordinates_transformation(sph_cphio, C_orbit)

# transform sph to cartesian (SIII)
domain_SIII = spherical_to_cartesian(sph_SIII)
# add time values in first column
domain_SIII = np.concatenate((time_array.reshape(-1, 1), domain_SIII), axis=1)



# ---------- Magnetic Field ---------

# jovian field
B_external = Bext_Community(domain_SIII.transpose())

# current sheet
B_sheet = B_sheet_Community(domain_SIII.transpose())
B_total = B_external + B_sheet

# induced field
r_core = 0.1 * R_C ; r_ocean = 0.9 * R_C ; r_surface = R_C ; r_iono = 1.05 * R_C
sig_core = 1e-9 ; sig_ocean = 10 ; sig_surface = 1e-9 ; sig_iono = 0.01
radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]

#B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, 3, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)
B_induced = B_induced_finite_conductivity_multilayer(domain_cphio.transpose(), B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

# total field
B_total = B_external + B_sheet + B_induced
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

