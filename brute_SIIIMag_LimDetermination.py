import matplotlib.pyplot as plt
import numpy as np
from functions import *

def get_general_data():
    '''
    Load spice data and return dictionary of 7d arrays storing values for t, x, y, z, r, theta, phi
    mission = "J" (juice) or "G" (galileo)
    '''
# define file names spice_data\callisto_wrt_jupiter_SIII_mag_longperiod1.csv
data_path = "./spice_data/callisto_wrt_jupiter_SIII_mag_longperiod1.csv"

# open spice data file
data = np.loadtxt(data_path, delimiter=",", unpack=True)
data = data.transpose()
Xjc, Yjc, Zjc, vxjc, vyjc, vzjc, t = data

# convert positions to m
cart = 1e3 * np.array([Xjc, Yjc, Zjc]).transpose()

# and to spherical coords
spher = cartesian_to_spherical(cart)

# combine t, cart, spher arrays
z = np.c_[t, cart]
z = np.transpose(np.c_[z, spher])
print(max(np.array(z[3])/R_J)) #  4.49 R_J
print(min(np.array(z[3])/R_J)) # -4.46 R_J
'''
plt.figure()
plt.plot(z[0], np.array(z[3])/R_J)
plt.show()
'''