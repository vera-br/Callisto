# list of planetary constants used
import numpy as np

# Jupiter
R_J = 71492e3 # m
J_spin_period = 10.1 * 3600 # s - synodic period at callisto
J_omega = 2 * np.pi / J_spin_period # Hz

# Callisto
R_C = 2410.3e3 # m
C_spin_period = 16.689 * 24 * 3600 # s
C_omega = 2 * np.pi / C_spin_period # Hz

# Sun
R_S = 696340e3 # m

# Astronomical Unit
AU = 15e10 # m