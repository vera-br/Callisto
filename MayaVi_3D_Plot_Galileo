import numpy as np
from mayavi import mlab
from functions import *

# load orbit data
Galileo, _B = get_pds_data()

# 3D plotting section

common_lim = 10

# Make background white.
mlab.figure(bgcolor=(1, 1, 1))  
mlab.view(azimuth=45, elevation=135)    

# Draws transparent pipe spanning the desired size for the axes because otherwise axes only stretch to span the last plotted thing
axis = [-common_lim, common_lim]
axis = mlab.plot3d(axis, axis, axis, opacity=0, line_width=0.01, tube_radius=0.1, color=(1,1,1))

# Draws the axes
axes = mlab.axes(color=(0, 0, 0), ranges=(-common_lim, common_lim, -common_lim, common_lim, -common_lim, common_lim), nb_labels=5)

axes.title_text_property.color = (0.0, 0.0, 0.0)
axes.title_text_property.font_family = 'times'

axes.label_text_property.color = (0.0, 0.0, 0.0)
axes.label_text_property.font_family = 'times'

axes.axes.font_factor = 1.0

axes.axes.label_format = '%-6.3g'

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
for orbit, vector in Galileo.items():
    x = vector[1] / R_C ; y = vector[2] / R_C ; z = vector[3] / R_C ; r = vector[4] / R_C
    
    # limiting coordinates to be inside the axes
    j = 0
    while abs(x[j]) > common_lim or abs(y[j]) > common_lim or abs(z[j]) > common_lim:
        j += 1

    k = j ; end_loop = False
    while abs(x[k]) < common_lim and abs(y[k]) < common_lim and abs(z[k]) < common_lim and end_loop == False:
        if k < len(r) - 1:
            k += 1
        else:
            end_loop = True

    x = x[j:k] ; y = y[j:k] ; z = z[j:k]

    # plotting the trajectories as tubes
    trajectory = mlab.plot3d(x, y, z,line_width=0.01, tube_radius=0.1, color=(colors[i][0], colors[i][1], colors[i][2]))
    i += 1

# makes size of objects independent from distance from the camera position
mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
# switches on axes orientation indicator
mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.
# shows plot
mlab.show()
