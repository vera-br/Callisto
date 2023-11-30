from functions import *


# load orbit data
Juice = get_spice_data('juice', 'callisto', 'cphio', 'J')
callisto_sun_jupsunorb = get_spice_data('callisto', 'sun', 'jupsunorb', 'J')
juice_callisto_jupsunorb = get_spice_data('juice', 'callisto', 'jupsunorb', 'J')
juice_cal_cphio_CA = get_closest_approach_data('juice', 'callisto', 'cphio', 'J')

common_lim = 1.1
create_callisto_plot(common_lim, ionosphere_CSO=True)
colors = colors_21()

times = {}
indices = {}
# plots all 21 orbits
for i in range(1,22): 
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
    r = vector[4] / R_C ; t = vector[0]
    spher_coords = cartesian_to_spherical(np.transpose([x,y,z])).transpose()
    theta = spher_coords[1] ; phi = spher_coords[2]
    
    indexs = [] ; ts = [] ; vec_x = [] ; vec_y = [] ; vec_z = []
    for h in range(len(y)):
        ionosphere_radius = 1.05 + (1 - 2 * np.arccos(np.sin(phi[h]) * np.sin(theta[h])) / np.pi) * 0.05
        if r[h] < ionosphere_radius:
                ts.append(t[h]) ; vec_x.append(x[h]) ; vec_y.append(y[h]) ; vec_z.append(z[h]) ; indexs.append(h)

    times['orbit%s' % (i)] = ts
    indices['orbit%s' % (i)] = indexs
    if len(vec_x) > 0:
        arrow_pos = [x[min_index], y[min_index], z[min_index]]
        arrow_vector = np.array([x[min_index+1]-x[min_index], y[min_index+1]-y[min_index], z[min_index+1]-z[min_index]])
        arrow_unit_vector = arrow_vector / np.linalg.norm(arrow_vector)

        # plotting the trajectories as tubes
        tube_radius = 0.025
        trajectory = mlab.plot3d(vec_x, vec_y, vec_z,line_width=0.01, tube_radius=tube_radius, color=(colors[i][0], colors[i][1], colors[i][2]))
        arrow = mlab.quiver3d(arrow_pos[0], arrow_pos[1], arrow_pos[2], arrow_unit_vector[0], arrow_unit_vector[1], arrow_unit_vector[2], line_width=0.1, color=(colors[i][0], colors[i][1], colors[i][2]), mode='cone')
        arrow_height = tube_radius * 10
        arrow.glyph.glyph_source.glyph_source.height = arrow_height
        arrow.glyph.glyph_source.glyph_source.center = np.array([arrow_height / 2, 0.  , 0.  ])
        arrow.glyph.glyph_source.glyph_source.radius = tube_radius * 2
# shows plot
mlab.show()

times2 = {}
for orbit, time in times.items():
    if len(time) > 0:
        t = max(time) - min(time)
        times2[str(orbit)] = t
print(times2)

