import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

orbit_n = 5

galileo_jupiter_SIII_spice = get_spice_data('juice', 'jupiter', 'SIII', 'J')
gal_jup_SIII_n_spice = galileo_jupiter_SIII_spice['orbit%s' % (orbit_n)]
t_spice = gal_jup_SIII_n_spice[0] / R_J
x_spice = gal_jup_SIII_n_spice[1] / R_J
y_spice = gal_jup_SIII_n_spice[2] / R_J
z_spice = gal_jup_SIII_n_spice[3] / R_J
r_spice = gal_jup_SIII_n_spice[4] / R_J

jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=False)
jm.Con2020.Config(equation_type='analytic', CartesianIn=True, CartesianOut=False)

Bx_spice, By_spice, Bz_spice = jm.Internal.Field(x_spice, y_spice, z_spice)
Bx_spiced, By_spiced, Bz_spiced = jm.Con2020.Field(x_spice, y_spice, z_spice)
B_mag_spice = np.sqrt(Bx_spice**2 + By_spice**2 + Bz_spice**2)
B_mag_spiced = np.sqrt(Bx_spiced**2 + By_spiced**2 + Bz_spiced**2)
Bx_spicet = Bx_spice + Bx_spiced
By_spicet = By_spice + By_spiced
Bz_spicet = Bz_spice + Bz_spiced
B_mag_spicet = np.sqrt(Bx_spicet**2 + By_spicet**2 + Bz_spicet**2)

fig, ax = plt.subplots(1,4)
ax[0].plot(t_spice, Bx_spicet, 'k')
ax[0].plot(t_spice, Bx_spiced, 'b')
ax[0].plot(t_spice, Bx_spice, 'r')

ax[1].plot(t_spice, By_spicet, 'k')
ax[1].plot(t_spice, By_spiced, 'b')
ax[1].plot(t_spice, By_spice, 'r')

ax[2].plot(t_spice, Bz_spice, 'r')
ax[2].plot(t_spice, Bz_spiced, 'b')
ax[2].plot(t_spice, Bz_spicet, 'k')

ax[3].plot(t_spice, B_mag_spicet, 'k', label='Total')
ax[3].plot(t_spice, B_mag_spiced, 'b', label='Disc')
ax[3].plot(t_spice, B_mag_spice, 'r', label='Internal')
ax[3].legend()

plt.show()