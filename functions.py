# File containing functions for data processing and plotting

# load modules
import numpy as np
import math as math
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as patches
import warnings
from mayavi import mlab
#import mayavi

# define constants
R_C = 2410.3 * 1e3
R_J = 71492 * 1e3
R_S = 696340 * 1e3
AU = 150000000 * 1e3

# coordinate systems

def cartesian_to_spherical(coords): 
    '''
    Changes 3d array of Cartesian coordinates (x,y,z) into Spherical coords (r,theta,phi)
    '''
    r = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    phi = np.sign(coords[:,1]) * np.arccos(coords[:,0] / (np.sqrt(coords[:,0]**2 + coords[:,1]**2)))
    theta = np.arccos(coords[:,2]/r)

    return np.array([r, theta, phi]).transpose()

def cartesian_to_spherical_single(coords): 
    '''
    Changes 3d array of Cartesian coordinates (x,y,z) into Spherical coords (r,theta,phi)
    '''
    r = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
    phi = np.sign(coords[1]) * np.arccos(coords[0] / (np.sqrt(coords[0]**2 + coords[1]**2)))
    theta = np.arccos(coords[2]/r)

    return np.array([r, theta, phi]).transpose()

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def cylindrical_coords(x, y, z):
    rho = np.sqrt(x**2 + y**2)
    phi = phi = np.sign(y) * np.arccos(x / (np.sqrt(x**2 + y**2)))
    return rho, phi, z

def CSunO_find_axis_unit_vectors(theta, phi):
    x_hat = np.array([-np.sin(phi), np.cos(phi), 0])
    y_hat = np.array([-np.sin(theta) * np.cos(phi), -np.sin(theta) * np.sin(phi), -np.cos(theta)])
    z_hat = np.array([-np.cos(theta) * np.cos(phi), -np.cos(theta) * np.sin(phi), np.sin(theta)])
    return (x_hat, y_hat, z_hat)

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def find_nearest_index(array, value):
    idx = np.searchsorted(array, value, side="left")
    return idx

def find_nearest_trajectories_G(target, reference_point, frame):
    '''
    finds trajectory from spice data with timestamps closest to those of the Galileo pds data
    '''
    Galileo, _ = get_pds_data()
    trajectories = get_spice_data(target, reference_point, frame, 'G')
    closest_trajectories = {}
    for i in range(0,7):
        G_vector = Galileo['orbit%s' % (i+1)]
        vector = trajectories['orbit%s' % (i+1)]
        trajectory = []
        for j in range(len(G_vector[0])):
            index = find_nearest_index(vector[0], G_vector[0][j])
            point = vector[:,int(index)]
            trajectory.append(point)
        closest_trajectories['orbit%s' % (i+1)] = np.transpose(trajectory)
    return closest_trajectories

# loading data from spice kernels

def get_spice_data(target, reference_point, frame, mission):
    '''
    Load spice data and return dictionary of 7d arrays storing values for t, x, y, z, r, theta, phi
    mission = "J" (juice) or "G" (galileo)
    '''
    # define file names
    data_path = "./spice_data/" + target + "_wrt_" + reference_point + "_" + frame + "_" + mission
    data_path_all = []

    #create array with all file names
    if mission == "J":
        for i in range(1, 22):
            data_path_all.append(data_path + str(i) + ".csv")
    else:
        for i in range(1, 8):
            data_path_all.append(data_path + str(i) + ".csv")

    # create dictionary for orbit data
    orbits_all = {}

    for i in range(len(data_path_all)):

        # open spice data file
        data = np.loadtxt(data_path_all[i], delimiter=",", unpack=True)
        data = data.transpose()
        Xjc, Yjc, Zjc, vxjc, vyjc, vzjc, t = data

        # convert positions to m
        cart = 1e3 * np.array([Xjc, Yjc, Zjc]).transpose()

        # and to spherical coords
        spher = cartesian_to_spherical(cart)

        # combine t, cart, spher arrays
        z = np.c_[t, cart]
        z = np.c_[z, spher]

        # save array as dictionary item
        orbits_all['orbit%s' % (i+1)] = np.transpose(z)

    return orbits_all


def get_pds_data():
    '''
    Load PDS data and return two dictionaries with 7d (t, x, y, z, r, theta, phi) and 5d arrays (t, Bx, By, Bz, |B|)
    '''
    orbits_all = {}
    bfield_all = {}

    data_path = "./galileo-mag-jup-calibrated/galileo_wrt_callisto_cphio_G"
    data_path_all = []

    #create array with all file names
    for i in range(1, 8):
        data_path_all.append(data_path + str(i) + ".TAB")

    for i in range(len(data_path_all)):

        t=[]

        # open data file
        data = np.loadtxt(data_path_all[i], usecols=(1, 2, 3, 4, 5, 6, 7))
        time_str = np.loadtxt(data_path_all[i], usecols=0, dtype=str)
        
        # Extract data components
        Bx = data[:, 0]
        By = data[:, 1]
        Bz = data[:, 2]
        B_mag = data[:, 3]
        x = data[:, 4] * R_C
        y = data[:, 5] * R_C
        z = data[:, 6] * R_C

        # combine arrays
        cart= np.c_[x, y]
        cart= np.c_[cart, z]

        # convert to spherical coords
        spher = cartesian_to_spherical(cart)

        # convert string to unix
        for x in range(len(time_str)):
            time = datetime.strptime(time_str[x], "%Y-%m-%dT%H:%M:%S.%f")        
            # Calculate the difference in seconds from a reference point
            timestamp_difference = time - datetime(2000, 1, 1, 12, 00)
            # Store the difference as a floating-point number
            t.append(timestamp_difference.total_seconds())

        # combine t, cart, spher arrays
        u = np.c_[t, cart]
        u = np.c_[u, spher]

        # combine t, bfield arrays
        B = np.c_[t, Bx]
        B = np.c_[B, By]
        B = np.c_[B, Bz]
        B = np.c_[B, B_mag]

        # save arrays in dictionaries
        orbits_all['orbit%s' % (i+1)] = np.transpose(u)
        bfield_all['bfield%s' % (i+1)] = np.transpose(B)

    return orbits_all, bfield_all

def closest_approach_data_G(target, reference_point, frame, mission):
    '''
    Compute CA data for a dictionary of orbits and return in its own dictionary of 7d arrays
    '''    
    CA_time_all =[]

    galileo_CA_data = {}
    CA_data_all = {}
    
    Galileo, Galileo_meas = get_pds_data()

    i = 0
    for key, array in Galileo.items():
        i += 1
        vector = np.transpose(array)
        min_index = np.argmin(vector[:, 4])
        CA_time_all.append(vector[min_index, 0])
        vector_CA = vector[min_index]
        vector_CA = np.append(vector_CA, min_index)
        galileo_CA_data['CA_orbit%s' % (i)] = vector_CA

    if target == "galileo":
        return galileo_CA_data
        
    else:
        dictionary = get_spice_data(target, reference_point, frame, mission)

        i = 0
        for key, array in dictionary.items():
            vector = np.transpose(array)
            time = vector[:, 0]
            CA = find_nearest(time, CA_time_all[i])
            index = int(np.where(time == CA)[0])
            i += 1
            CA_data_all['CA_orbit%s' % (i)] = vector[index]        
            
        return CA_data_all

def CA_info(orbit):
    '''
    input: orbit as a vector array
    returns: CA vector (t, x, y, z, r, theta, phi, min_index)
    '''
    # finds index of smallest r value
    min_index = np.argmin(orbit[4])

    # finds full vector at min. index
    CA_info_vector_i = np.transpose(orbit)[min_index,:] 

    # appends min. index to CA vector for use in other functions
    CA_info_vector_i = np.append(CA_info_vector_i, min_index)

    return CA_info_vector_i

juice_cal_cphio_CA = 0

gal_cal_cphio_CA = 0

def get_closest_approach_data(target, reference_point, frame, mission):
    '''
    returns dictionary of closest approaches
    '''
    global juice_cal_cphio_CA
    global gal_cal_cphio_CA
    
    if mission == 'J': 
        # if juice_callisto_cphio closest approach not already calculated
        if juice_cal_cphio_CA == False:
            juice_cal_cphio_CA = {}

            # gets full orbit info. for juice_callisto_cphio
            orbits_all_jcalcphio = get_spice_data('juice', 'callisto', 'cphio', 'J')

            i = 1
            for orbit, vector in orbits_all_jcalcphio.items():
                # finds closest approaches for juice_callisto_cphio and adds to dictionary
                CA_info_vector = CA_info(vector)
                juice_cal_cphio_CA['CA_orbit%s' %(i)] = CA_info_vector
                i += 1

        if target == 'juice' and reference_point == 'callisto' and frame == 'cphio':
            return juice_cal_cphio_CA
        
        orbits_all_i = get_spice_data(target, reference_point, frame, mission)
        closest_approach_vectors_i = {}

        i = 1
        for orbit, vector in orbits_all_i.items():
            orbit_i = np.transpose(vector)

            # takes min. index included in dictionary of juice_callisto_cphio CA vectors
            minindex = int(juice_cal_cphio_CA['CA_orbit%s' % (i)][7])

            # finds closest approach vector for respective bodies and frame with min. index
            closest_approach_vectors_i['CA_orbit%s' % (i)] = orbit_i[minindex]
            i += 1

    
    if mission == 'G':
            
        # if galileo_callisto_cphio closest approach not already calculated
        if gal_cal_cphio_CA == False:
            gal_cal_cphio_CA = {}

            # gets full orbit info. for galileo_callisto_cphio
            orbits_all_galcalcphio, _gal_B = get_pds_data()

            i = 1
            for orbit, vector in orbits_all_galcalcphio.items():
                # finds closest approaches for galileo_callisto_cphio and adds to dictionary
                CA_info_vector = CA_info(vector)
                gal_cal_cphio_CA['CA_orbit%s' %(i)] = CA_info_vector
                i += 1

        if target == 'galileo':
            return gal_cal_cphio_CA
        
        orbits_all_i = get_spice_data(target, reference_point, frame, mission)
        closest_approach_vectors_i = {}

        i = 1
        for orbit, vector in orbits_all_i.items():
            orbit_i = np.transpose(vector)

            # takes min. time included in dictionary of galileo_callisto_cphio CA vectors
            mintime = gal_cal_cphio_CA['CA_orbit%s' % (i)][0]

            # finds closest time from spice data to get vector at closest approach
            vector = np.transpose(vector)
            time = vector[:, 0]
            CA = find_nearest(time, mintime)
            index = int(np.where(time == CA)[0])
            closest_approach_vectors_i['CA_orbit%s' % (i)] = vector[index]

            i += 1

    return closest_approach_vectors_i

def create_callisto_plot(common_lim, ionosphere_CSO=False, night_cone_CSO=False):
    # Make background white.
    scene = mlab.figure(bgcolor=(1, 1, 1))  

    # Draws transparent pipe spanning the desired size for the axes because otherwise axes only stretch to span the last plotted thing
    axis = [-common_lim, common_lim]
    axis = mlab.plot3d(axis, axis, axis, opacity=0, line_width=0.01, tube_radius=0.0001, color=(1,1,1))

    # Draws the axes
    axes = mlab.axes(color=(0, 0, 0), ranges=(-common_lim, common_lim, -common_lim, common_lim, -common_lim, common_lim), nb_labels=5)

    axes.title_text_property.color = (0.0, 0.0, 0.0)
    axes.title_text_property.font_family = 'times'

    axes.label_text_property.color = (0.0, 0.0, 0.0)
    axes.label_text_property.font_family = 'times'

    axes.axes.font_factor = 1.0

    axes.axes.label_format = '%-6.3g'

    mlab.outline(color=(0, 0, 0))

    # makes size of objects independent from distance from the camera position
    mlab.gcf().scene.parallel_projection = True  # Source: <<https://stackoverflow.com/a/32531283/2729627>>.
    # switches on axes orientation indicator
    mlab.orientation_axes()  # Source: <<https://stackoverflow.com/a/26036154/2729627>>.

    plot_callisto()

    if ionosphere_CSO == True:
        plot_ionosphere_CSO()
    
    if night_cone_CSO == True:
        night_cone = mlab.quiver3d(0,0,0, 0, -1, 0, color=(0,0,0), mode='cylinder')
        night_cone.glyph.glyph_source.glyph_source.radius = 1.0
        night_cone.glyph.glyph_source.glyph_source.resolution = 256
        night_cone.glyph.glyph_source.glyph_source.height = common_lim
        night_cone.glyph.glyph_source.glyph_source.center = np.array([ 0.  , -0.5 * common_lim,  0.  ])
        night_cone.actor.property.edge_tint = np.array([1., 1., 1.])
        night_cone.actor.property.emissive_factor = np.array([1., 1., 1.])
        night_cone.actor.property.selection_color = np.array([1., 0., 0., 1.])
        night_cone.actor.property.opacity = 0.1
        


def plot_callisto():
    sphere = mlab.points3d(0,0,0, color=(0,0,0), resolution=256)
    sphere.glyph.glyph_source.glyph_source.radius = 1

def plot_ionosphere_CSO():
    phis = np.linspace(0, 2 * np.pi, 181)
    thetas = np.linspace(0, np.pi, 91)

    phis_grid, thetas_grid = np.meshgrid(phis, thetas)
    x_grid = np.zeros_like(phis_grid)
    y_grid = np.zeros_like(phis_grid)
    z_grid = np.zeros_like(phis_grid)

    for i in range(len(thetas)):
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

    ionosphere = mlab.mesh(x_grid, y_grid, z_grid, color=(0,0,1), opacity=0.2)

def colors_7():
    # defines 7 equally spaced colours around edge of colour wheel
    colors_x = [1, 1, 2/7, 0, 0, 2/7, 1]
    colors_y = [0, 6/7, 1, 1, 4/7, 0, 0]
    colors_z = [0, 0, 0, 4/7, 1, 1, 6/7]
    colors = np.array([colors_x, colors_y, colors_z])
    colors = np.transpose(colors)
    return colors

def colors_21():
    # defines 21 equally spaced colours around edge of colour wheel
    colors_x = [1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1]
    colors_y = [0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7, 0, 0, 0, 0, 0, 0, 0]
    colors_z = [0, 0, 0, 0, 0, 0, 0, 0, 2/7, 4/7, 6/7, 1, 1, 1, 1, 1, 1, 1, 6/7, 4/7, 2/7]
    colors = np.array([colors_x, colors_y, colors_z])
    colors = np.transpose(colors)
    return colors