from trajectories.trajectory_analysis import *

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')
juice_cal_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

# 3D plotting section

common_lim = 10
create_callisto_plot(common_lim)
colors = colors_21()
# plots all 21 orbits

for i in range(1,22): 
    # dayside group = C4-9 requires range(4,10), nightside group = C13-17 requires range(13,18)
    vector = Juice['orbit%s'%(i)]
    juice_cal_cphio_CA_i = juice_cal_cphio_CA['CA_orbit%s' % (i)]
    min_index = int(juice_cal_cphio_CA_i[7])

    x = vector[1] / R_C ; y = vector[2] / R_C ; z = vector[3] / R_C ; r = vector[4] / R_C
    arrow_pos = [x[min_index], y[min_index], z[min_index]]
    arrow_vector = np.array([x[min_index+1]-x[min_index], y[min_index+1]-y[min_index], z[min_index+1]-z[min_index]])
    
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
    arrow_unit_vector = arrow_vector / np.linalg.norm(arrow_vector)

    tube_radius = 0.1
    trajectory = mlab.plot3d(x, y, z, line_width=0.01, tube_radius=tube_radius, color=(colors[i-1][0], colors[i-1][1], colors[i-1][2]))
    arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=0.1, color=(colors[i-1][0], colors[i-1][1], colors[i-1][2]), mode='cone')
    arrow_height = tube_radius * 10
    arrow.glyph.glyph_source.glyph_source.height = arrow_height
    arrow.glyph.glyph_source.glyph_source.center = np.array([arrow_height / 2, 0.  , 0.  ])
    arrow.glyph.glyph_source.glyph_source.radius = tube_radius * 2
    print(i)
    i += 1

# shows plot
mlab.show()
