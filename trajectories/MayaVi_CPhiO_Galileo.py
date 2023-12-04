from trajectories.trajectory_analysis import *

# load orbit data
Galileo, _B = get_pds_data()
galileo_cal_cphio_CA  = get_closest_approach_data('galileo', 'callisto', 'cphio', 'G')

# 3D plotting section

common_lim = 10

create_callisto_plot(common_lim)
colors = colors_7()

# plots all 21 orbits
i = 0
for orbit, vector in Galileo.items():
    galileo_cal_cphio_CA_i = galileo_cal_cphio_CA['CA_orbit%s' % (i+1)]
    min_index = int(galileo_cal_cphio_CA_i[7])
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
    arrow_pos = [x[min_index], y[min_index], z[min_index]]
    arrow_vector = np.array([x[min_index+1]-x[min_index], y[min_index+1]-y[min_index], z[min_index+1]-z[min_index]])
    arrow_unit_vector = arrow_vector / np.linalg.norm(arrow_vector)

    tube_radius = 0.1
    trajectory = mlab.plot3d(x, y, z, line_width=0.01, tube_radius=tube_radius, color=(colors[i][0], colors[i][1], colors[i][2]))
    arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=0.1, color=(colors[i][0], colors[i][1], colors[i][2]), mode='cone')
    arrow_height = tube_radius * 10
    arrow.glyph.glyph_source.glyph_source.height = arrow_height
    arrow.glyph.glyph_source.glyph_source.center = np.array([arrow_height / 2, 0.  , 0.  ])
    arrow.glyph.glyph_source.glyph_source.radius = tube_radius * 2

    i += 1

# shows plot
mlab.show()
