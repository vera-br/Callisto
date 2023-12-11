import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

_, B_gal = get_pds_data()
t_G = _['orbit1'][0]
x_G = _['orbit1'][1] * 100 / R_J
y_G = _['orbit1'][2] * 100 / R_J
z_G = _['orbit1'][3] * 100 / R_J

B_gal_1 = B_gal['bfield1']

galileo_jupiter_SIII_GK = get_spice_data('galileo', 'jupiter', 'SIII', 'GK')
galileo_jup_SIII_1 = galileo_jupiter_SIII_GK['orbit1']
t = galileo_jup_SIII_1[0]
x = galileo_jup_SIII_1[1] / R_J
y = galileo_jup_SIII_1[2] / R_J
z = galileo_jup_SIII_1[3] / R_J

Bx, By, Bz = jm.Internal.Field(x, y, z)
B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
Bx_G = B_gal_1[1] ; By_G = B_gal_1[2] ; Bz_G = B_gal_1[3]
B_mag_G = np.sqrt(Bx_G**2 + By_G**2 + Bz_G**2)

print('min t_G = ' + str(min(t_G)))
print('max t_G = ' + str(max(t_G)))
print('min t = ' + str(min(t)))
print('max t = ' + str(max(t)))

fig, ax = plt.subplots(1, 4)

ax[0].plot(t, Bx)

ax[1].plot(t, By)

ax[2].plot(t, Bz)

ax[3].plot(t, B_mag)

ax[0].plot(t_G, Bx_G)

ax[1].plot(t_G, By_G)

ax[2].plot(t_G, Bz_G)

ax[3].plot(t_G, B_mag_G)

plt.show()

plt.figure()

plt.plot(t, Bx)

plt.plot(t, By)

plt.plot(t, Bz)

plt.plot(t, B_mag)

plt.show()

fig, ax = plt.subplots(1, 3)

ax[0].plot(t, x)

ax[1].plot(t, y)

ax[2].plot(t, z)

ax[0].plot(t_G, x_G)

ax[1].plot(t_G, y_G)

ax[2].plot(t_G, z_G)


plt.show()

