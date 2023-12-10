# functions to load data from spice and pds, as well as cartesial to spherical conversions

import numpy as np
from datetime import datetime

# constants
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

        # convert times to timestamps
        #t = [Timestamp(datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')) for timestamp in t]

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


# load data from pds

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
            '''
            time = Timestamp(time_str[x])
            t.append(time)
            '''
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