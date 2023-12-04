from functions import *

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_sun_jupsunorb = get_spice_data('callisto', 'sun', 'jupsunorb', 'J')
juice_callisto_jupsunorb = get_spice_data('juice', 'callisto', 'jupsunorb', 'J')
juice_cal_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

# 3D plotting section

common_lim = 2
create_callisto_plot(common_lim, ionosphere_CSO=True)
colors = colors_21()

# plots all 21 orbits
i = 0

for i in range(1, 22): 
    # dayside group = C4-9 requires range(4,10), nightside group = C13-17 requires range(13,18)
    vector = juice_callisto_jupsunorb['orbit%s' % (i)]
    calsun_i = callisto_sun_jupsunorb['orbit%s' % (i)]
    juice_cal_cphio_CA_i = juice_cal_cphio_CA['CA_orbit%s' % (i)]
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
        tube_radius = 0.05
        trajectory = mlab.plot3d(x, y, z,line_width=0.01, tube_radius=0.05, color=(colors[i][0], colors[i][1], colors[i][2]))
        if min_index > j and min_index < k:
            arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=0.1, color=(colors[i][0], colors[i][1], colors[i][2]), mode='cone')
            arrow_height = tube_radius * 10
            arrow.glyph.glyph_source.glyph_source.height = arrow_height
            arrow.glyph.glyph_source.glyph_source.center = np.array([arrow_height / 2, 0.  , 0.  ])
            arrow.glyph.glyph_source.glyph_source.radius = tube_radius * 2
        
    i += 1

# shows plot
mlab.show()
