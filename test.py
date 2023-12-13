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
B_external = Bext_full(orbit_SIII)

# current sheet
rho, z, B_sheet = B_disk(orbit_SIII_mag, R0=R0, R1=R1, D=D, I_constant=Icon, azimuthal_field=True)

# induced field
B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit4_cphio, B_external, 0.03, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)

# total field
B_total = B_external + B_sheet + B_induced

#---------plot-----------
plot_time_evolution(B_total, orbit_cphio, orbit_CA, flyby_n, "Total")