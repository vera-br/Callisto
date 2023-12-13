import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

orbit_n = 1
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

gal_jup_SIII_calc = Galileo_trajectories_SIII_from_CPhiO()
gal_jup_SIII_n = gal_jup_SIII_calc['orbit%s' % (orbit_n)]
t_calc = gal_jup_SIII_n[0] / R_J
x_calc = gal_jup_SIII_n[1] / R_J
y_calc = gal_jup_SIII_n[2] / R_J
z_calc = gal_jup_SIII_n[3] / R_J
r_calc = gal_jup_SIII_n[4] / R_J

jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=True)

Bx_calc, By_calc, Bz_calc = jm.Internal.Field(x_calc, y_calc, z_calc)
B_mag_calc = np.sqrt(Bx_calc**2 + By_calc**2 + Bz_calc**2)

galileo_jupiter_SIII_spice = get_spice_data('galileo', 'jupiter', 'SIII', 'GK')
gal_jup_SIII_n_spice = galileo_jupiter_SIII_spice['orbit%s' % (orbit_n)]
t_spice = gal_jup_SIII_n_spice[0] / R_J
x_spice = gal_jup_SIII_n_spice[1] / R_J
y_spice = gal_jup_SIII_n_spice[2] / R_J
z_spice = gal_jup_SIII_n_spice[3] / R_J
r_spice = gal_jup_SIII_n_spice[4] / R_J

Bx_spice, By_spice, Bz_spice = jm.Internal.Field(x_spice, y_spice, z_spice)
B_mag_spice = np.sqrt(Bx_spice**2 + By_spice**2 + Bz_spice**2)

fig, ax = plt.subplots(1, 4)
ax[0].plot(t_calc, x_calc)
ax[0].plot(t_spice, x_spice)
ax[1].plot(t_calc, y_calc)
ax[1].plot(t_spice, y_spice)
ax[2].plot(t_calc, z_calc)
ax[2].plot(t_spice, z_spice)
ax[3].plot(t_calc, r_calc)
ax[3].plot(t_spice, r_spice)
plt.show()

fig, ax = plt.subplots(1, 4)

ax[0].plot(t_calc, Bx_G)
ax[0].plot(t_calc, Bx_calc)
ax[0].plot(t_spice, Bx_spice)

ax[1].plot(t_calc, By_G)
ax[1].plot(t_calc, By_calc)
ax[1].plot(t_spice, By_spice)

ax[2].plot(t_calc, Bz_G)
ax[2].plot(t_calc, Bz_calc)
ax[2].plot(t_spice, Bz_spice)

ax[3].plot(t_calc, B_mag_G)
ax[3].plot(t_calc, B_mag_calc)
ax[3].plot(t_spice, B_mag_spice)

plt.show()

galileo_callisto_cphio_spice = get_spice_data('galileo', 'callisto', 'cphio', 'GK')
gal_cal_cphio_spice = galileo_callisto_cphio_spice['orbit%s' % (orbit_n)]
t_cs = gal_cal_cphio_spice[0]
x_cs = gal_cal_cphio_spice[1] / R_C
y_cs = gal_cal_cphio_spice[2] / R_C
z_cs = gal_cal_cphio_spice[3] / R_C
r_cs = gal_cal_cphio_spice[4] / R_C

fig, ax = plt.subplots(1, 4)

ax[0].plot(t_cs, x_cs)
ax[0].plot(t_G, x_G)
ax[1].plot(t_cs, y_cs)
ax[1].plot(t_G, y_G)
ax[2].plot(t_cs, z_cs)
ax[2].plot(t_G, z_G)
ax[3].plot(t_cs, r_cs)
ax[3].plot(t_G, r_G)

plt.show()



# fig, ax = plt.subplots(1, 4)

# ax[0].plot(t_spice, Bx_spice)

# ax[1].plot(t_spice, By_spice)

# ax[2].plot(t_spice, Bz_spice)

# ax[3].plot(t_spice, B_mag_spice)

# plt.show()

# plt.figure()

# plt.plot(t_spice, Bx_spice)

# plt.plot(t_spice, By_spice)

# plt.plot(t_spice, Bz_spice)

# plt.plot(t_spice, B_mag_spice)

# plt.show()


