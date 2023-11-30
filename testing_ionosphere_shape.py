import numpy as np
import matplotlib.pyplot as plt
from mayavi.api import Engine
import mayavi
from mayavi import mlab
from functions import *
import pandas as pd

phis = np.linspace(0, 2 * np.pi, 181)
thetas = np.linspace(0, np.pi, 91)

phis_grid, thetas_grid = np.meshgrid(phis, thetas)
x_grid = np.zeros_like(phis_grid)
y_grid = np.zeros_like(phis_grid)
z_grid = np.zeros_like(phis_grid)

for i in range(len(thetas)):
    print(i)
    for j in range(len(phis)):
        h = j

        phi = phis_grid[i,j]
        theta = thetas_grid[i,j]
        r = 1.05 + (1 - 2 * np.arccos(np.sin(phi) * np.sin(theta)) / np.pi) * 0.05
        x = r * np.cos(phi) * np.sin(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(theta)
        x_grid[i,j] = x
        y_grid[i,j] = y
        z_grid[i,j] = z
        

common_lim = 1.1

# Make background white.
scene = mlab.figure(bgcolor=(1, 1, 1))  

# Draws transparent pipe spanning the desired size for the axes because otherwise axes only stretch to span the last plotted thing
axis = [-common_lim, common_lim]
axis = mlab.plot3d(axis, axis, axis, opacity=0, line_width=0.01, tube_radius=0.00001, color=(1,1,1))

# Draws the axes
axes = mlab.axes(color=(0, 0, 0), ranges=(-common_lim, common_lim, -common_lim, common_lim, -common_lim, common_lim), nb_labels=13)

axes.title_text_property.color = (0.0, 0.0, 0.0)
axes.title_text_property.font_family = 'times'

axes.label_text_property.color = (0.0, 0.0, 0.0)
axes.label_text_property.font_family = 'times'

axes.axes.font_factor = 1.0

axes.axes.label_format = '%-6.3g'

mlab.outline(color=(0, 0, 0))

radius = 1
sphere = mlab.points3d(0,0,0, color=(0,0,0), resolution=256, scale_factor=2*radius)
ionosphere = mlab.mesh(x_grid, y_grid, z_grid, color=(0,0,1), opacity=0.2)

# makes size of objects independent from distance from the camera position
mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
# switches on axes orientation indicator
mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.
# shows plot
mlab.show()
