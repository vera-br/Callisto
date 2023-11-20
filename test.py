from mayavi import mlab
from functions import *

sun_wrt_callisto_cphio_CA = get_closest_approach_data("sun", "callisto", "cphio", "J")

# get sun-cal angle
J_sun_cphio_vector = []
for orbit, vector in sun_wrt_callisto_cphio_CA.items():
    J_sun_cphio_vector.append(vector[1:4])
J_sun_cphio_vector = np.array(J_sun_cphio_vector)

# arrays for plotting day/night boundary
mlab.clf()
xx, yy = np.mgrid[-10:10:20j, -10:10:20j]

#plotting day/night boundary plane
zz = (-J_sun_cphio_vector[4, 0] * xx - J_sun_cphio_vector[4, 1] * yy) * 1. /J_sun_cphio_vector[4, 2]
plane = mlab.surf(zz, warp_scale='auto', color=(.5, .5, 0.5))

# makes size of objects independent from distance from the camera position
mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
# switches on axes orientation indicator
mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.
# shows plot
mlab.show()