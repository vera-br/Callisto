import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

_, B_gal = get_pds_data()
t_G = _['orbit1'][0]
x_G = _['orbit1'][1] / R_C
y_G = _['orbit1'][2] / R_C
z_G = _['orbit1'][3] / R_C

B_gal_1 = B_gal['bfield1']

galileo_jupiter_SIII_GK = get_spice_data('galileo', 'jupiter', 'SIII', 'GK')
galileo_jup_SIII_1 = galileo_jupiter_SIII_GK['orbit1']
t = galileo_jup_SIII_1[0]
x = galileo_jup_SIII_1[1] / R_J
y = galileo_jup_SIII_1[2] / R_J
z = galileo_jup_SIII_1[3] / R_J

galileo_callisto_cphio_GK = get_spice_data('galileo', 'callisto', 'cphio', 'GK')
gal_cal_cphio_1 = galileo_callisto_cphio_GK['orbit1']
t_c = gal_cal_cphio_1[0]
x_c = gal_cal_cphio_1[1] / R_C
y_c = gal_cal_cphio_1[2] / R_C
z_c = gal_cal_cphio_1[3] / R_C

galileo_callisto_IAUCAL_GK = get_spice_data('galileo', 'callisto', 'IAU_CALLISTO', 'GK')
gal_cal_IAUCAL_1 = galileo_callisto_IAUCAL_GK['orbit1']
t_i = gal_cal_IAUCAL_1[0]
x_i = gal_cal_IAUCAL_1[1] / R_C
y_i = gal_cal_IAUCAL_1[2] / R_C
z_i = gal_cal_IAUCAL_1[3] / R_C

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

ax[0].plot(t_c, x_c)

ax[1].plot(t_c, y_c)

ax[2].plot(t_c, z_c)

ax[0].plot(t_G, x_G)

ax[1].plot(t_G, y_G)

ax[2].plot(t_G, z_G)

ax[0].plot(t_i, x_i)

ax[1].plot(t_i, y_i)

ax[2].plot(t_i, z_i)


plt.show()

plt.figure()

r_c = np.sqrt(x_c**2 + y_c**2 + z_c**2)
r_i = np.sqrt(x_i**2 + y_i**2 + z_i**2)
r_G = _['orbit1'][4] / R_C
plt.plot(t_c, r_c)
plt.plot(t_G, r_G)
plt.plot(t_i, r_i)
plt.show()

