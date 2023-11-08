import numpy as np
from mayavi import mlab
from functions import *

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio')
Callisto = get_spice_data("callisto", "jupiter", "jupsunorb")
Callisto_SIII = get_spice_data("callisto", "jupiter", "SIII")
Jupiter = get_spice_data('jupiter', 'sun', 'IAU_SUN')

Galileo, Galileo_meas = get_pds_data()

# save CA data
Juice_CA, Callisto_CA, Callisto_SIII_CA, Jupiter_CA = closest_approach_data_4(Juice, Callisto, Callisto_SIII, Jupiter)
Galileo_CA = closest_approach_data(Galileo)

common_lim = 10

mlab.figure(bgcolor=(1, 1, 1))  # Make background white.

# Draws transparent pipe spanning the desired size for the axes because otherwise axes only stretch to span the next plotted thing
axis = [-common_lim, common_lim]
axis = mlab.plot3d(axis, axis, axis, opacity=0, line_width=0.01, tube_radius=0.1, color=(1,1,1))

# Draws the axes
axes = mlab.axes(color=(0, 0, 0), ranges=(-common_lim, common_lim, -common_lim, common_lim, -common_lim, common_lim), nb_labels=5)
axes.title_text_property.color = (0.0, 0.0, 0.0)
axes.title_text_property.font_family = 'times'
axes.label_text_property.color = (0.0, 0.0, 0.0)
axes.label_text_property.font_family = 'times'
mlab.outline(color=(0, 0, 0))

# plots sphere of specified radius
radius = 1
sphere = mlab.points3d(0,0,0, color=(0,0,0), resolution=256, scale_factor=2*radius)

# defines 21 equally spaced colours around edge of colour wheel
colors_x = [1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1]
colors_y = [0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0]
colors_z = [0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7]
colors = np.array([colors_x, colors_y, colors_z])
colors = np.transpose(colors)

# plots all 21 orbits
i = 0
for orbit, vector in Juice.items():
    x = vector[1] / R_C ; y = vector[2] / R_C ; z = vector[3] / R_C 
    trajectory = mlab.plot3d(x, y, z,line_width=0.01, tube_radius=0.1, color=(colors[i][0], colors[i][1], colors[i][2]))
    i += 1

# makes size of objects independent from distance from the camera position
mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
# switches on axes orientation indicator
mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.
# shows plot
mlab.show()
