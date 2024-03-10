# load modules and functions
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
from pandas import Timestamp
import numpy.fft as fft

from trajectories.trajectory_analysis import *
from plot_scripts.plot_field_components import *
from induced_field import *
from jupiter_field import *
from current_sheet import *
from field_functions import *
from khurana1997 import *

# load data

callisto_jupiter_SIII = get_spice_data_longperiod('callisto', 'jupiter', 'SIII', 'G')
callisto_jupiter_SIII_mag = get_spice_data_longperiod('callisto', 'jupiter', 'SIII_mag', 'G')
callisto_jupiter_JSO = get_spice_data_longperiod('callisto', 'jupiter', 'jupsunorb', 'G')


# specify orbit
flyby_n = 2

orbit_cal_SIII = callisto_jupiter_SIII["orbit%s" % (flyby_n)]
orbit_cal_SIII_mag = callisto_jupiter_SIII_mag["orbit%s" % (flyby_n)]
orbit_cal_JSO = callisto_jupiter_JSO["orbit%s" % (flyby_n)]

#---------magnetic fields-----------

# SIII coords
x0 = -33.5 #* R_J
rho0 = 33.2 #* R_J
v0 = 37.4 #* R_J hr^-1
C1 = 80.3 ; C2 = 690.4 ; C3 = 101.3 ; C4 = -1.7
a1 = 2.49 ; a2 = 1.80 ; a3 = 2.64
r01 = 38.0 #* R_J
rho02 = 2.14 #* R_J
rho03 = 12.5 #* R_J 
D1 = 2.01 #* R_J 
D2 = 13.27 #* R_J
p = 6.26e-3 ; q = 0.35
const_ = [x0, rho0, v0, C1, C2, C3, C4, a1, a2, a3, r01, rho02, rho03, D1, D2, p, q]
const_scale = [1.49801966, 1.23520174, 1.49650863, 1.43940198, 0.54528406, 1.31496966,
 1.493787,   1.45328466, 1.35288508, 0.88258155, 0.63155589, 0.73154669,
 0.8343032,  0.63042957, 0.54541552, 0.50317376, 0.53259648,]
consts = np.multiply(const_scale, const_)

B_external_cal = Bext_Community(orbit_cal_SIII)
B_sheet_cal = B_sheet_khurana(orbit_cal_JSO, orbit_cal_SIII_mag, orbit_cal_SIII)

B_full_ext_cal = B_external_cal + B_sheet_cal


fig, ax = plt.subplots(3,1)
ax[0].plot(orbit_cal_SIII[6] * 180/np.pi, B_full_ext_cal[:,0], 'or')
ax[1].plot(orbit_cal_SIII[6] * 180/np.pi, B_full_ext_cal[:,1], 'or')
ax[2].plot(orbit_cal_SIII[6] * 180/np.pi, B_full_ext_cal[:,2], 'or')

plt.show()

ts = np.mean([orbit_cal_SIII[0][i+1] - orbit_cal_SIII[0][i] for i in range(len(orbit_cal_SIII[0]) - 1)])
sr = 1 / ts
t = orbit_cal_SIII[0]


# Write a function DFT(x) which takes in one argument, 
# x - input 1 dimensional real-valued signal. 
# The function will calculate the DFT of the signal and return the DFT values. 

def DFT(x):
    """
    Function to calculate the 
    discrete Fourier Transform 
    of a 1D real-valued signal x
    """

    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * np.pi * k * n / N)
    
    X = np.dot(e, x)
    
    return X

# X = DFT(np.transpose(B_full_ext_cal)[0])
X = fft.fft(np.transpose(B_full_ext_cal)[0])

# calculate the frequency
N = len(X)
n = np.arange(N)
T = N/sr
freq = n/T 

n_oneside = N//2
# get the one side frequency
f_oneside = freq[:n_oneside]

# normalize the amplitude
X_oneside =X[:n_oneside]/n_oneside

fig, ax = plt.subplots()
ax.stem(1/f_oneside / 3600, abs(X_oneside), 'r', markerfmt=' ')
ax.set_xscale('log')
# ax.set_yscale('log')
plt.show()

