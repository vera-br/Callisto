import numpy as np
import matplotlib.pyplot as plt
from mayavi.api import Engine
import mayavi
from mayavi import mlab
from functions import *
import pandas as pd

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_jupiter_jupsunorb = get_spice_data('callisto', 'jupiter', 'jupsunorb', 'J')

callisto_sun_jupsunorb = get_spice_data('callisto', 'sun', 'jupsunorb', 'J')
juice_callisto_jupsunorb = get_spice_data('juice', 'callisto', 'jupsunorb', 'J')

juice_cal_cphio_CA = {}

i = 1
for orbit, vector in Juice.items():
    # finds closest approaches for juice_callisto_cphio and adds to dictionary
    CA_info_vector = CA_info(vector)
    juice_cal_cphio_CA['CA_orbit%s' %(i)] = CA_info_vector
    i += 1

def CSunO_find_axis_unit_vectors(theta, phi):
    x_hat = np.array([-np.sin(phi), np.cos(phi), 0])
    y_hat = np.array([-np.sin(theta) * np.cos(phi), -np.sin(theta) * np.sin(phi), -np.cos(theta)])
    z_hat = np.array([-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)])
    return (x_hat, y_hat, z_hat)

# 3D plotting section

common_lim = 1.5

# Make background white.
scene = mlab.figure(bgcolor=(1, 1, 1))  

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
ionosphere = mlab.points3d(0,0,0, resolution=256, opacity=0.2, color=(0,0,1))
ionosphere.glyph.glyph_source.glyph_source.radius = 1.2

# defines 21 equally spaced colours around edge of colour wheel
colors_x = [1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1]
colors_y = [0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0]
colors_z = [0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7]
colors = np.array([colors_x, colors_y, colors_z])
colors = np.transpose(colors)

# plots all 21 orbits
i = 0

juice_cal_CSO_CA = {}

for orbit, vector in juice_callisto_jupsunorb.items():
    calsun_i = callisto_sun_jupsunorb['orbit%s' % (i+1)]
    juice_cal_cphio_CA_i = juice_cal_cphio_CA['CA_orbit%s' % (i+1)]
    min_index = int(juice_cal_cphio_CA_i[7])
    thetas = calsun_i[5]
    phis = calsun_i[6]
    x_new = [] ; y_new = [] ; z_new = []
    for m in range(len(vector[1])):
        x_hat, y_hat, z_hat = CSunO_find_axis_unit_vectors(thetas[m], phis[m])
        new_x = np.dot(x_hat, np.transpose(vector[1:4, :])[m]) ; x_new.append(new_x)
        new_y = np.dot(y_hat, np.transpose(vector[1:4, :])[m]) ; y_new.append(new_y)
        new_z = np.dot(z_hat, np.transpose(vector[1:4, :])[m]) ; z_new.append(new_z)

    x = np.array(x_new) / R_C ; y = np.array(y_new) / R_C ; z = np.array(z_new) / R_C
    
    t_CA = juice_cal_cphio_CA_i[0]
    x_CA = np.dot(x_hat, np.transpose(vector[1:4, :])[min_index])
    y_CA = np.dot(y_hat, np.transpose(vector[1:4, :])[min_index]) 
    z_CA = np.dot(z_hat, np.transpose(vector[1:4, :])[min_index])
    spher_coords_CA = cartesian_to_spherical_single([x_CA, y_CA, z_CA])
    CA_array = [t_CA, x_CA, y_CA, z_CA]
    CA_array.append(spher_coords_CA[0])
    CA_array.append(spher_coords_CA[1])
    CA_array.append(spher_coords_CA[2])
    juice_cal_CSO_CA['CA_orbit%s' % (i+1)] = CA_array

    # limiting coordinates to be inside the axes
    j = 0 ; j_out_of_bounds = False
    while abs(x[j]) > common_lim or abs(y[j]) > common_lim or abs(z[j]) > common_lim:
        if j < len(x) - 1:
            j += 1
        else:
            j_out_of_bounds = True
            break

    k = j ; end_loop2 = False
    while abs(x[k]) < common_lim and abs(y[k]) < common_lim and abs(z[k]) < common_lim and j_out_of_bounds == False and end_loop2 == False:
        if k < len(x) - 1:
            k += 1
        else:
            end_loop2 = True

    if j_out_of_bounds == False:
        arrow_pos = [x[min_index], y[min_index], z[min_index]]
        arrow_vector = np.array([x[min_index+1]-x[min_index], y[min_index+1]-y[min_index], z[min_index+1]-z[min_index]])
        arrow_unit_vector = arrow_vector / np.linalg.norm(arrow_vector)

        x = x[j:k] ; y = y[j:k] ; z = z[j:k]

        # plotting the trajectories as tubes
        trajectory = mlab.plot3d(x, y, z,line_width=0.01, tube_radius=0.1, color=(colors[i][0], colors[i][1], colors[i][2]))
        arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=2, color=(colors[i][0], colors[i][1], colors[i][2]), mode='cone')
    i += 1

pd.DataFrame(juice_cal_CSO_CA).to_csv('juice_cal_CSO.csv', index=False)

night_cone = mlab.quiver3d(0,0,0, 0, -1, 0, color=(0,0,0), mode='cylinder')
night_cone.glyph.glyph_source.glyph_source.radius = 1.0
night_cone.glyph.glyph_source.glyph_source.resolution = 256
night_cone.glyph.glyph_source.glyph_source.height = 1.5
night_cone.glyph.glyph_source.glyph_source.center = np.array([ 0.  , -0.75,  0.  ])
night_cone.actor.property.edge_tint = np.array([1., 1., 1.])
night_cone.actor.property.emissive_factor = np.array([1., 1., 1.])
night_cone.actor.property.selection_color = np.array([1., 0., 0., 1.])
night_cone.actor.property.opacity = 0.1

# makes size of objects independent from distance from the camera position
mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
# switches on axes orientation indicator
mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.
# shows plot
mlab.show()
