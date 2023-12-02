"""
The following script calculates the magnetic field components for a given trajectory
"""
from jupiter_external import *
from induced_field import *


# constants
pi = constants.pi
mu0 = constants.mu_0
RJ = 71492e3  # Jupiter radius
RG = 2634.1e3  # Ganymede radius

### Load in data from SPICE kernels, data loaded in has units km and km/s
# data paths
juice_ganymede_data_path = ""
juice_jupiter_data_path = ""
ganymede_jupiter_data_path = ""

# juice to ganymede
juice_wrt_ganymede = np.loadtxt(
    juice_ganymede_data_path, delimiter=",", unpack=True
)
juice_wrt_ganymede = juice_wrt_ganymede.transpose()
Xjg, Yjg, Zjg, vxjg, vyjg, vzjg, t = juice_wrt_ganymede

# juice to jupiter
juice_wrt_jupiter = np.loadtxt(
    juice_jupiter_data_path, delimiter=",", unpack=True
)
juice_wrt_jupiter = juice_wrt_jupiter.transpose()
Xjj, Yjj, Zjj, vxjj, vyjj, vzjj, _ = juice_wrt_jupiter

# ganymede to jupiter
ganymede_wrt_jupiter = np.loadtxt(
    ganymede_jupiter_data_path, delimiter=",", unpack=True
)
ganymede_wrt_jupiter = ganymede_wrt_jupiter.transpose()
Xgj, Ygj, Zgj, vxgj, vygj, vzgj, _ = ganymede_wrt_jupiter

# convert positions to m
juice_wrt_ganymede = 1e3 * np.array([Xjg, Yjg, Zjg]).transpose()
juice_wrt_jupiter = 1e3 * np.array([Xjj, Yjj, Zjj]).transpose()
ganymede_wrt_jupiter = 1e3 * np.array([Xgj, Ygj, Zgj]).transpose()

# set initial time to zero
t = t - t[0]

# induced field parameters
rm = 2634.1e3  # Ganymede radius
r0 = rm  # outer conducting layer radius
r1 = rm - 200e3  # inner conducting layer radius
sigma = 0.1

# rotation speeds
J_spin_period = 9.9 * 3600  # actual period is 9.9 hours
G_orbit_period = 172 * 3600
G_omega = 2 * pi / G_orbit_period  # hours^-1
J_omega = 2 * pi / J_spin_period  # hours^-1

# induced bessel parameters
omega_synodic = 2 * pi / J_spin_period
omega_eccentric = 2 * pi / G_orbit_period

# parameters for cylindrical plasma sheet
R0 = 7.8  # disc inner radius (RJ)
R1 = 51.4  # disc outer radius (RJ)
D = 3.6  # disc half thickness (RJ)
Icon = 139.6  # current constant = mu0 * I / 2 (nT)
thetaD = np.radians(9.3)  # disc normal from rotation axis (radians)
phiD = np.radians(204.2)  # azimuth angle of disc normal (radians)

# calculate magnetict fields
B_intrinsic = B_intrinsic_vectors(positions=juice_wrt_ganymede)
B_external = Bext_full(positions=juice_wrt_jupiter, times=t)

Bext_ecc = Bext_eccentric_variation(positions=ganymede_wrt_jupiter)
Bind_ecc_superconductor = B_induced_superconductor(
    pos_vectors=juice_wrt_ganymede, Bext_vectors=Bext_ecc, r0=r0, rm=rm
)
Bind_ecc_finite = B_induced_finite_conductivity(
    pos_vectors=juice_wrt_ganymede,
    Bext_vectors=Bext_ecc,
    sigma=0.05,
    omega=omega_synodic,
    rm=rm,
    r0=r0,
    r1=r1,
)

Bext_syn = Bext_synodic_variation(pos_avg=np.mean(juice_wrt_jupiter, axis=0), times=t)
Bind_syn_superconductor = B_induced_superconductor(
    pos_vectors=juice_wrt_ganymede, Bext_vectors=Bext_syn, r0=r0, rm=rm
)
Bind_syn_finite = B_induced_finite_conductivity(
    pos_vectors=juice_wrt_ganymede,
    Bext_vectors=Bext_syn,
    sigma=0.05,
    omega=omega_synodic,
    rm=rm,
    r0=r0,
    r1=r1,
)

z1, rho1, B_plasma_sheet1 = B_sheet(
    positions=ganymede_wrt_jupiter / RJ,
    times=t,
    R0=R0,
    R1=R1,
    D=D,
    I_constant=Icon,
    thetaD=thetaD,
)
z2, rho2, B_plasma_sheet2 = B_sheet(
    positions=juice_wrt_jupiter / RJ,
    times=t,
    R0=R0,
    R1=R1,
    D=D,
    I_constant=Icon,
    thetaD=thetaD,
)