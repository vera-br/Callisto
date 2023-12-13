import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

orbit_n = 1

# PDS data for Galileo's position and B field in CPhiO
vectors_gal, B_gal = get_pds_data()
B_gal_n = B_gal['bfield%s' % (orbit_n)]
vector_gal_n = vectors_gal['orbit%s' % (orbit_n)]
t_G = vector_gal_n[0]
x_G = vector_gal_n[1] / R_C
y_G = vector_gal_n[2] / R_C
z_G = vector_gal_n[3] / R_C
r_G = vector_gal_n[4] / R_C
Bx_G = B_gal_n[1] ; By_G = B_gal_n[2] ; Bz_G = B_gal_n[3]
B_mag_G = np.sqrt(Bx_G**2 + By_G**2 + Bz_G**2)

# spice data for Galileo's position in CPhiO
galileo_callisto_cphio_spice = get_spice_data('galileo', 'callisto', 'cphio', 'GK')
gal_cal_cphio_spice = galileo_callisto_cphio_spice['orbit%s' % (orbit_n)]
t_cs = gal_cal_cphio_spice[0]
x_cs = gal_cal_cphio_spice[1] / R_C
y_cs = gal_cal_cphio_spice[2] / R_C
z_cs = gal_cal_cphio_spice[3] / R_C
r_cs = gal_cal_cphio_spice[4] / R_C

# initialises magnetic field model used to calculate B field
# set to take Cartesian input and output in Spherical coordinates
jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=False)
jm.Con2020.Config(equation_type='analytic', CartesianIn=True, CartesianOut=False)

# calculates Galileo's position in SIII
gal_jup_SIII_calc = Galileo_trajectories_SIII_from_CPhiO()
gal_jup_SIII_n = gal_jup_SIII_calc['orbit%s' % (orbit_n)]
t_calc = gal_jup_SIII_n[0] / R_J
x_calc = gal_jup_SIII_n[1] / R_J
y_calc = gal_jup_SIII_n[2] / R_J
z_calc = gal_jup_SIII_n[3] / R_J
r_calc = gal_jup_SIII_n[4] / R_J

# calculates Jupiter's B field from Galileo's calculated SIII position
Br_calc, Btheta_calc, Bphi_calc = jm.Internal.Field(x_calc, y_calc, z_calc)
Br_calcd, Btheta_calcd, Bphi_calcd = jm.Con2020.Field(x_calc, y_calc, z_calc)
B_mag_calc = np.sqrt(Br_calc**2 + Btheta_calc**2 + Bphi_calc**2)
B_mag_calcd = np.sqrt(Br_calcd**2 + Btheta_calcd**2 + Bphi_calcd**2)
Br_calct = Br_calc + Br_calcd
Btheta_calct = Btheta_calc + Btheta_calcd
Bphi_calct = Bphi_calc + Bphi_calcd
B_mag_calct = np.sqrt(Br_calct**2 + Btheta_calct**2 + Bphi_calct**2)

# finds Galileo's position in SIII from saved spice kernel data
galileo_jupiter_SIII_spice = get_spice_data('galileo', 'jupiter', 'SIII', 'GK')
gal_jup_SIII_n_spice = galileo_jupiter_SIII_spice['orbit%s' % (orbit_n)]
t_spice = gal_jup_SIII_n_spice[0] / R_J
x_spice = gal_jup_SIII_n_spice[1] / R_J
y_spice = gal_jup_SIII_n_spice[2] / R_J
z_spice = gal_jup_SIII_n_spice[3] / R_J
r_spice = gal_jup_SIII_n_spice[4] / R_J

# calculates Jupiter's B field from Galileo's position from spice kernel
Br_spice, Btheta_spice, Bphi_spice = jm.Internal.Field(x_spice, y_spice, z_spice)
Br_spiced, Btheta_spiced, Bphi_spiced = jm.Con2020.Field(x_spice, y_spice, z_spice)
B_mag_spice = np.sqrt(Br_spice**2 + Btheta_spice**2 + Bphi_spice**2)
B_mag_spiced = np.sqrt(Br_spiced**2 + Btheta_spiced**2 + Bphi_spiced**2)
Br_spicet = Br_spice + Br_spiced
Btheta_spicet = Btheta_spice + Btheta_spiced
Bphi_spicet = Bphi_spice + Bphi_spiced
B_mag_spicet = np.sqrt(Br_spicet**2 + Btheta_spicet**2 + Bphi_spicet**2)

# converts from spherical SIII to cartesian CPhiO
Bx_calc = Bphi_calc ; By_calc = -Br_calc ; Bz_calc = -Btheta_calc
Bx_calcd = Bphi_calcd ; By_calcd = -Br_calcd ; Bz_calcd = -Btheta_calcd
Bx_calct = Bphi_calct ; By_calct = -Br_calct ; Bz_calct = -Btheta_calct
Bx_spice = Bphi_spice ; By_spice = -Br_spice ; Bz_spice = -Btheta_spice
Bx_spiced = Bphi_spiced ; By_spiced = -Br_spiced ; Bz_spiced = -Btheta_spiced
Bx_spicet = Bphi_spicet ; By_spicet = -Br_spicet ; Bz_spicet = -Btheta_spicet

# plots Jupiter's internal (red), magnetodisc (blue) and total (black) magnetic field at Galileo
fig, ax = plt.subplots(1,4)

fig.suptitle('')

ax[0].set_title('B_x')
ax[0].plot(t_spice, Bx_spice, 'r')
ax[0].plot(t_spice, Bx_spiced, 'b')
ax[0].plot(t_spice, Bx_spicet, 'k')

ax[1].set_title('B_y')
ax[1].plot(t_spice, By_spice, 'r')
ax[1].plot(t_spice, By_spiced, 'b')
ax[1].plot(t_spice, By_spicet, 'k')

ax[2].set_title('B_z')
ax[2].plot(t_spice, Bz_spice, 'r')
ax[2].plot(t_spice, Bz_spiced, 'b')
ax[2].plot(t_spice, Bz_spicet, 'k')

ax[3].set_title('|B|')
ax[3].plot(t_spice, B_mag_spice, 'r', label='Internal')
ax[3].plot(t_spice, B_mag_spiced, 'b', label='Disc')
ax[3].plot(t_spice, B_mag_spicet, 'k', label='Total')
ax[3].legend()

plt.show()

# compares Galileo's calculated (red) to spice (blue) position
fig, ax = plt.subplots(1, 4)

ax[0].set_title('x')
ax[0].plot(t_spice, x_spice, 'b')
ax[0].plot(t_calc, x_calc, 'r')


ax[1].set_title('y')
ax[1].plot(t_spice, y_spice, 'b')
ax[1].plot(t_calc, y_calc, 'r')


ax[2].set_title('z')
ax[2].plot(t_spice, z_spice, 'b')
ax[2].plot(t_calc, z_calc, 'r')


ax[3].set_title('r')
ax[3].plot(t_spice, r_spice, 'b', label='Spice')
ax[3].plot(t_calc, r_calc, 'r', label='Calculated')
ax[3].legend()
plt.show()

# compares PDS B field to total Jupiter field at calculated SIII position
fig, ax = plt.subplots(1, 4)

ax[0].set_title('B_x')
ax[0].plot(t_calc, Bx_G, 'r')
ax[0].plot(t_calc, Bx_calct, 'b')

ax[1].set_title('B_y')
ax[1].plot(t_calc, By_G, 'r')
ax[1].plot(t_calc, By_calct, 'b')

ax[2].set_title('B_z')
ax[2].plot(t_calc, Bz_G, 'r')
ax[2].plot(t_calc, Bz_calct, 'b')

ax[3].set_title('|B|')
ax[3].plot(t_calc, B_mag_G, 'r', label='PDS')
ax[3].plot(t_calc, B_mag_calct, 'b', label='Calculated')
ax[3].legend()

plt.show()

# comparing Galileo's spice (red) position to PDS (blue)

fig, ax = plt.subplots(1, 4)

ax[0].set_title('x')
ax[0].plot(t_cs, x_cs, 'r')
ax[0].plot(t_G, x_G, 'b')

ax[1].set_title('y')
ax[1].plot(t_cs, y_cs, 'r')
ax[1].plot(t_G, y_G, 'b')

ax[2].set_title('z')
ax[2].plot(t_cs, z_cs, 'r')
ax[2].plot(t_G, z_G, 'b')

ax[3].set_title('r')
ax[3].plot(t_cs, r_cs, 'r', label='Spice')
ax[3].plot(t_G, r_G, 'b', label='PDS')
ax[3].legend()

plt.show()



