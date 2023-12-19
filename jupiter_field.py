# jovian field

from sph_harmonics.coefficients import *
import sph_harmonics.legendre as legendre
from scipy.spatial.transform import Rotation as ROT
import JupiterMag as jm

# matthew and ciaran's files (that we want to get rid of?)
from maths import unit_spherical_to_cartesian
from sph_harmonics.magf_higher_order import *


# constants
R_J = 71492e3
J_spin_period = 10.1 * 3600
C_spin_period = 16.689 * 24 * 3600
C_omega = 2 * np.pi / C_spin_period
J_omega = 2 * np.pi / J_spin_period

# Spherical harmonic coefficients
quad_coeffs = Jupiter_quad_coeffs_2022()

# limits of summation
m = np.arange(0, 3, 1)
n = np.arange(1, 3, 1)

def Bext_full(orbit):
    """
    Calculates the external field vectors due to the eccentric orbit and synodic rotation of Jupiter.
    :param orbit: array with t (J200), x (m), y(m), z(m), r(m), theta(deg), phi(deg) 
    :return: time evolution array of Bx, By, Bz in nT
    """
    # time intervals and angles
    rot_theta = J_omega * orbit[0]

    orbit = orbit.transpose()

    # calculate Jupiter's field at Ganymede
    magfJ_vectors = []
    for angle, vector in zip(rot_theta, orbit):

        pos = vector[1:4]       
        r = vector[4]
        theta = vector[5]
        phi = vector[6]

        # Calculate legendre polynomials over a range of theta and create dictionary
        leg, dleg = legendre.quadrupole(theta)

        # calculate magnetic field components of Jupiter
        Br, Btheta, Bphi = B_sph_components(
            r=r,
            a=R_J,
            theta=theta,
            phi=phi,
            coeffs=quad_coeffs,
            n_ints=n,
            m_ints=m,
            leg_poly=leg,
            dleg_poly=dleg,
        )

        # Convert spherical to cartesian coordinates
        Bsph = np.array([Br, Btheta, Bphi])
        Bcart = unit_spherical_to_cartesian(Bsph, theta, phi)
        Bx, By, Bz = Bcart

        rotation = ROT.from_rotvec(angle * np.array([0, 0, 1]))

        # arrange point coordinates in shape (N, 3) for vectorized processing
        B_XYZ = np.array([Bx.ravel(), By.ravel(), Bz.ravel()]).transpose()

        # apply rotation
        B_XYZrot = rotation.apply(B_XYZ)

        # return to original shape of meshgrid
        BXrot = B_XYZrot[:, 0].reshape(Bx.shape)
        BYrot = B_XYZrot[:, 1].reshape(By.shape)
        BZrot = B_XYZrot[:, 2].reshape(Bz.shape)
        B_cart_rot = np.array([BXrot, BYrot, BZrot])

        magfJ_vectors.append(B_cart_rot)
    return np.array(magfJ_vectors)

def Bext_Community(orbit_SIII):
    jm.Internal.Config(Model='jrm33', CartesianIn=True, CartesianOut=False)
    x = orbit_SIII[1] / R_J
    y = orbit_SIII[2] / R_J
    z = orbit_SIII[3] / R_J
    Br, Btheta, Bphi = jm.Internal.Field(x, y, z)
    Bx = Bphi
    By = -Br
    Bz = -Btheta
    return np.array([Bx, By, Bz]).transpose()
