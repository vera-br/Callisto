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
from field_functions import *

# load data
galileo_wrt_callisto_cphio, B_PDSs = get_pds_data()
galileo_wrt_callisto_cphio_CA = get_closest_approach_data("galileo", "callisto", "cphio", "G")
galileo_wrt_jupiter_SIII = Galileo_trajectories_SIII_from_CPhiO()



# specify orbit
flyby_n = 2

orbit_cphio = galileo_wrt_callisto_cphio["orbit%s" % (flyby_n)]
orbit_SIII = galileo_wrt_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_CA = galileo_wrt_callisto_cphio_CA["CA_orbit%s" % (flyby_n)]
B_PDS = B_PDSs['bfield%s' % (flyby_n)]
B_mag = np.sqrt(B_PDS[1]**2 + B_PDS[2]**2 + B_PDS[3]**2)

#---------magnetic fields-----------

# jovian field
# B_external = Bext_full(orbit_SIII)
B_external = Bext_Community(orbit_SIII)


# current sheet
# rho, z, B_sheet = B_disk(orbit_SIII_mag, R0=R0, R1=R1, D=D, I_constant=Icon, azimuthal_field=True)
B_sheet = B_sheet_Community(orbit_SIII)
B_total = B_external + B_sheet

# induced field
r_core = 0.1 * R_C ; r_ocean = 0.9 * R_C ; r_surface = R_C ; r_iono = 1.05 * R_C
sig_core = 1e-9 ; sig_ocean = 10 ; sig_surface = 1e-9 ; sig_iono = 1e-9
radii = [r_core, r_ocean, r_surface, r_iono]
conductivities = [sig_core, sig_ocean, sig_surface, sig_iono]

#B_induced = B_induced_infinite(orbit_cphio, B_external + B_sheet, R_C, R_C - 100)
#B_induced = B_induced_finite_conductivity(orbit_cphio, B_external, 3, 2*np.pi /(10.1*3600), R_C, R_C - 100, R_C - 200)
B_induced = B_induced_finite_conductivity_multilayer(orbit_cphio, B_external + B_sheet, 2*np.pi /(10.1*3600), conductivities, radii)

# total field
B_total = B_external + B_sheet + B_induced
B_mag_tot = np.sqrt(B_total[:, 0]**2 + B_total[:, 1]**2 + B_total[:, 2]**2)

#---------plot-----------
plot_time_evolution_Gal(B_total, orbit_cphio, orbit_CA, flyby_n, "Total")

fig, ax = plt.subplots(2, 2)
ax[0,0].plot(B_PDS[0], B_PDS[1], label='PDS')
ax[0,0].plot(orbit_cphio[0], B_total[:, 0], label='Calc.')
ax[0,0].set_title('Bx')

ax[0,1].plot(B_PDS[0], B_PDS[1], label='PDS')
ax[0,1].plot(orbit_cphio[0], B_total[:, 1], label='Calc.')
ax[0,1].set_title('By')

ax[1,0].plot(B_PDS[0], B_PDS[1], label='PDS')
ax[1,0].plot(orbit_cphio[0], B_total[:, 2], label='Calc.')
ax[1,0].set_title('Bz')

ax[1,1].plot(B_PDS[0], B_mag, label='PDS')
ax[1,1].plot(orbit_cphio[0], B_mag_tot, label='Calc.')
ax[1,1].set_title('|B|')
ax[1,1].legend()

fig.suptitle('r_ocean = ' + str(r_ocean / R_C) + ', r_iono = ' + str(r_iono / R_C) + ', sig_ocean = ' + str(sig_ocean) + ', sig_iono = ' + str(sig_iono))
plt.show()
