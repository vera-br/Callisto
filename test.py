# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *

# load data
juice_wrt_callisto_cphio = get_spice_data("juice", "callisto", "cphio", "J")
juice_wrt_callisto_cphio_CA = get_closest_approach_data("juice", "callisto", "cphio", "J")
juice_wrt_jupiter_SIII = get_spice_data("juice", "jupiter", "SIII", "J")
juice_wrt_jupiter_SIII_mag = get_spice_data("juice", "jupiter", "SIII_mag", "J")
#callisto_wrt_jupiter_SIII = get_spice_data("callisto", "jupiter", "SIII", "J")


# specify orbit
flyby_n = 5

orbit_cphio = juice_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = juice_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_SIII_mag = juice_wrt_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_CA = juice_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]


#---------magnetic fields-----------

# jovian field
# B_external = Bext_full(orbit_SIII)
B_external = Bext_Community(orbit_SIII)


# current sheet
# rho, z, B_sheet = B_disk(orbit_SIII_mag, R0=R0, R1=R1, D=D, I_constant=Icon, azimuthal_field=True)
B_sheet = B_sheet_Community(orbit_SIII)
B_total = B_external + B_sheet

# induced field
r_core = 0.1 * R_C ; r_ocean = 0.9 * R_C ; r_surface = R_C ; r_iono = 1.1 * R_C
sig_core = 1e-9 ; sig_ocean = 10 ; sig_surface = 1e-9 ; sig_iono = 0.1
radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]
#B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, 3, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)
B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

# total field
B_total = B_external + B_sheet + B_induced

#---------plot-----------
plot_time_evolution(B_induced, orbit_cphio, orbit_CA, flyby_n, "Induced")