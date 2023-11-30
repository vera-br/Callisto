import numpy as np
from mayavi import mlab
from functions import *

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')

# 3D plotting section

common_lim = 10

create_callisto_plot(common_lim)
colors = colors_21()

# plots all 21 orbits
i = 0
for orbit, vector in Juice.items():
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

# shows plot
mlab.show()
