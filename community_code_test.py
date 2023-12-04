import JupiterMag as jm
from trajectories.trajectory_analysis import *
import matplotlib.pyplot as plt

juice_jup_SIII = get_spice_data('juice', 'jupiter', 'SIII', 'J')
juice_jup_SIII_5 = juice_jup_SIII['orbit5']
t = juice_jup_SIII_5[0]
x = juice_jup_SIII_5[1]
y = juice_jup_SIII_5[2]
z = juice_jup_SIII_5[3]

jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=True)

Bx, By, Bz = jm.Internal.Field(x, y, z)
B_mag = abs(Bx) + abs(By) + abs(Bz)

fig, ax = plt.subplots(1, 4)

ax[0].plot(t, Bx)
ax[1].plot(t, By)
ax[2].plot(t, Bz)
ax[3].plot(t, B_mag)

plt.show()

jm.Internal.Test()
plt.show()


