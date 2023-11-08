# File containing functions for data processing and plotting

# load modules
import numpy as np
from datetime import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as patches
import warnings

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
    phi = np.sign(coords[:,0]) * np.arccos(coords[:,1] / (np.sqrt(coords[:,0]**2 + coords[:,1]**2)))
    # theta = np.nan_to_num(x=theta)    - can use this for nan error if necessary
    theta = np.arccos(coords[:,2]/r)

    return np.array([r, theta, phi]).transpose()


# loading data from spice kernels

def get_spice_data(target, reference_point, frame):
    '''
    Load spice data and return dictionary of 7d arrays storing values for t, x, y, z, r, theta, phi
    '''
    # define file names
    data_path = "./spice_data/" + target + "_wrt_" + reference_point + "_" + frame + "_J"
    data_path_all = []

    #create array with all file names
    for i in range(1, 22):
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

        # convert string to utc
        for i in range(len(time_str)):
            time = datetime.fromisoformat(time_str[i]).timestamp()
            t.append(time)

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

def closest_approach_data(dictionary):
    '''
    Compute CA data for a dictionary of orbits and return in its own dictionary of 7d arrays
    '''    
    CA_data = {}
    CA_data_all = {}
    i = 0

    for key, array in dictionary.items():
        i += 1
        vector = np.transpose(array)
        min_index = np.argmin(vector[:, 4])
        CA_data_all['CA_orbit%s' % (i)] = vector[min_index]
        
    return CA_data_all

def closest_approach_data_4(dictionary, dict2, dict3, dict4):
    '''
    Same as closest_approach_data but for 5 dictionaries
    '''    
    CA_data = {}
    CA_data_all = {}
    i = 0

    indices =[]
    CA_data_all_2 = {}
    CA_data_all_3 = {}
    CA_data_all_4 = {}
    CA_data_all_5 = {}

    for key, array in dictionary.items():
        i += 1
        vector = np.transpose(array)
        min_index = np.argmin(vector[:, 4])
        indices.append(min_index)
        CA_data_all['CA_orbit%s' % (i)] = vector[min_index]

    i = 0
    for key, array in dict2.items():
        
        vector = np.transpose(array)
        CA_data_all_2['CA_orbit%s' % (i)] = vector[indices[i]]
        i += 1

    i = 0
    for key, array in dict3.items():
        
        vector = np.transpose(array)
        CA_data_all_3['CA_orbit%s' % (i)] = vector[indices[i]]
        i += 1

    i = 0
    for key, array in dict4.items():
        
        vector = np.transpose(array)
        CA_data_all_4['CA_orbit%s' % (i)] = vector[indices[i]]
        i += 1
           
    return CA_data_all, CA_data_all_2, CA_data_all_3, CA_data_all_4

def CA_info(orbit):
    '''
    takes orbit as a vector
    '''
    # finds index of smallest r value
    min_index = np.argmin(orbit[4])

    # finds full vector at min. index
    CA_info_vector_i = np.transpose(orbit)[min_index,:] 

    # appends min. index to CA vector for use in other functions
    CA_info_vector_i = np.append(CA_info_vector_i, min_index)

    return CA_info_vector_i

def get_orbit_data(target, reference_point, frame):
    '''
    returns a dictionary of all orbits
    '''
    
    datapath = "./spice_data/" + target + "_wrt_" + reference_point + "_" + frame + "_J"
    data_paths_i = []

    for i in range(1, 22):
        data_paths_i.append(datapath + str(i) + ".csv")

    orbits_all_i = {}

    for i in range(len(data_paths_i)):
        t, cart, spher = data_processing(data_paths_i[i])
        z = np.c_[t, cart]
        z = np.c_[z, spher]
        orbits_all_i['orbit%s' % (i+1)] = np.transpose(z)
    return orbits_all_i

juice_cal_cphio_CA = 0

def get_closest_approach_data(target, reference_point, frame):
    '''
    returns dictionary of closest approaches
    '''
    global juice_cal_cphio_CA

    orbits_all_i = get_spice_data(target, reference_point, frame)
    closest_approach_vectors_i = {}
    
    # if juice_callisto_cphio closest approach not already calculated
    if juice_cal_cphio_CA == False:
        juice_cal_cphio_CA = {}

        # gets full orbit info. for juice_callisto_cphio
        orbits_all_jcalcphio = get_spice_data('juice', 'callisto', 'cphio')

        i = 1
        for orbit, vector in orbits_all_jcalcphio.items():
            # finds closest approaches for juice_callisto_cphio and adds to dictionary
            CA_info_vector = CA_info(vector)
            juice_cal_cphio_CA['CA_orbit%s' %(i)] = CA_info_vector
            i += 1

    i = 1
    for orbit, vector in orbits_all_i.items():
        orbit_i = np.transpose(vector)

        # takes min. index included in dictionary of juice_callisto_cphio CA vectors
        minindex = int(juice_cal_cphio_CA['CA_orbit%s' % (i)][7])

        # finds closest approach vector for respective bodies and frame with min. index
        closest_approach_vectors_i['CA_orbit%s' % (i)] = orbit_i[minindex]
        i += 1

    return closest_approach_vectors_i
# plots

#def plot_trajectories():
