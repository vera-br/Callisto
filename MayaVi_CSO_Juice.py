from functions import *

# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_sun_jupsunorb = get_spice_data('callisto', 'sun', 'jupsunorb', 'J')
juice_callisto_jupsunorb = get_spice_data('juice', 'callisto', 'jupsunorb', 'J')

juice_cal_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

# 3D plotting section

common_lim = 1.1

create_callisto_plot(common_lim, ionosphere_CSO=True, night_cone_CSO=True)
colors = colors_21()

# plots all 21 orbits
i = 0

juice_cal_CSO_CA = {}

for orbit, vector in juice_callisto_jupsunorb.items():
    calsun_i = callisto_sun_jupsunorb['orbit%s' % (i+1)]
    juice_cal_cphio_CA_i = juice_cal_cphio_CA['CA_orbit%s' % (i+1)]
    min_index = int(juice_cal_cphio_CA_i[7])
    thetas = calsun_i[5]
    print(len(thetas))
    phis = calsun_i[6]
    x_new = [] ; y_new = [] ; z_new = []
    for m in range(len(vector[1])):
        x_hat, y_hat, z_hat = CSunO_find_axis_unit_vectors(thetas[m], phis[m])
        new_x = np.dot(x_hat, np.transpose(vector[1:4, :])[m]) ; x_new.append(new_x)
        new_y = np.dot(y_hat, np.transpose(vector[1:4, :])[m]) ; y_new.append(new_y)
        new_z = np.dot(z_hat, np.transpose(vector[1:4, :])[m]) ; z_new.append(new_z)

    x = np.array(x_new) / R_C ; y = np.array(y_new) / R_C ; z = np.array(z_new) / R_C
    '''
    # collects CSO CA vector to save into pd array
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
    '''
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
        tube_radius = 0.025
        trajectory = mlab.plot3d(x, y, z,line_width=0.01, tube_radius=tube_radius, tube_sides=12, color=(colors[i][0], colors[i][1], colors[i][2]))
        if min_index > j and min_index < k:
            arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=0.1, color=(colors[i][0], colors[i][1], colors[i][2]), mode='cone')
            arrow_height = tube_radius * 10
            arrow.glyph.glyph_source.glyph_source.height = arrow_height
            arrow.glyph.glyph_source.glyph_source.center = np.array([arrow_height / 2, 0.  , 0.  ])
            arrow.glyph.glyph_source.glyph_source.radius = tube_radius * 2
    i += 1

# saves the CSO CA vector to .csv
# pd.DataFrame(juice_cal_CSO_CA).to_csv('juice_cal_CSO.csv', index=False)


# shows plot
mlab.show()
