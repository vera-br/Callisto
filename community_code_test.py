import JupiterMag as jm
from trajectories.trajectory_analysis import *
from field_functions import *
import matplotlib.pyplot as plt

_, B_gal = get_pds_data()
gal_jup_SIII = Galileo_trajectories_SIII_from_CPhiO()
juice_jup_SIII_5 = gal_jup_SIII['orbit1']
B_gal_1 = B_gal['bfield1']
t = juice_jup_SIII_5[0] / R_J
x = juice_jup_SIII_5[1] / R_J
y = juice_jup_SIII_5[2] / R_J
z = juice_jup_SIII_5[3] / R_J

jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=True)

Bx, By, Bz = jm.Internal.Field(x, y, z)
B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)
Bx_G = B_gal_1[1] ; By_G = B_gal_1[2] ; Bz_G = B_gal_1[3]
B_mag_G = np.sqrt(Bx_G**2 + By_G**2 + Bz_G**2)

fig, ax = plt.subplots(1, 4)

ax[0].plot(t, Bx)
ax[0].plot(t, B_gal_1[1])
ax[1].plot(t, By)
ax[1].plot(t, B_gal_1[2])
ax[2].plot(t, Bz)
ax[2].plot(t, B_gal_1[3])
ax[3].plot(t, B_mag)
ax[3].plot(t, B_mag_G)

plt.show()

jm.Internal.Test()
plt.show()


