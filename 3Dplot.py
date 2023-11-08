from functions import *
from itertools import cycle

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio')
Callisto = get_spice_data("callisto", "jupiter", "jupsunorb")
Callisto_SIII = get_spice_data("callisto", "jupiter", "SIII")
Jupiter = get_spice_data('jupiter', 'sun', 'IAU_SUN')

Galileo, Galileo_meas = get_pds_data()

# save CA data
Juice_CA, Callisto_CA, Callisto_SIII_CA, Jupiter_CA = closest_approach_data_4(Juice, Callisto, Callisto_SIII, Jupiter)
Galileo_CA = closest_approach_data(Galileo)

# Create a color palette with 21 different colors
colors = plt.cm.tab20.colors
color_cycle = cycle(colors)

# 3D Plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Defining and applying the common limits
lim = 3

common_xlim = (-lim, lim)  
common_ylim = (-lim, lim) 
common_zlim = (-lim, lim) 

ax.set_xlim(common_xlim)
ax.set_ylim(common_ylim)
ax.set_zlim(common_zlim)

# Set axis labels
ax.set_xlabel('x [$R_C$]')
ax.set_ylabel('y [$R_C$]')
ax.set_zlabel('z [$R_C$]')

# Make Callisto surface
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the Callisto surface
ax.plot_surface(x, y, z, color='grey')

# Plot the trajectories (all of them with loop)
for orbit, vector in Juice.items():
    color = next(color_cycle)

    ax.plot(vector[1] / R_C, vector[2] / R_C, vector[3] / R_C, color=color)

plt.show()