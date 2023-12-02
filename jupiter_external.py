import numpy as np
from scipy import constants
from maths import cartesian_to_spherical
from maths import unit_spherical_to_cartesian
from sph_harmonics.magf_higher_order import *
from sph_harmonics.coefficients import *
import sph_harmonics.legendre as legendre
from scipy.spatial.transform import Rotation as ROT

# constants
RJ = 71492e3  # Jupiter radius

# limits of summation
m = np.arange(0, 3, 1)
n = np.arange(1, 3, 1)

J_spin_period = 10.1 * 3600
C_spin_period = 17 * 24 * 3600
C_omega = 2 * np.pi / C_spin_period
J_omega = 2 * np.pi / J_spin_period
omega = J_omega - C_omega

# Spherical harmonic coefficients
quad_coeffs = Jupiter_quad_coeffs_2022()


def Bext_synodic_variation(times, pos_avg):
    """
    Calculates the external field vectors with variation due Jupiter's synodic rotation
    :param times:
    :param pos_avg: average position vector of callisto

    :return: array of external field vectors
    """
    # time intervals and angles
    rot_theta = omega * times

    # convert rotated cartesian mesh to spherical coordinates
    r, theta, phi = cartesian_to_spherical(pos_avg)

    # Calculate legendre polynomials over a range of theta and create dictionary
    leg, dleg = legendre.quadrupole(theta)

    # calculate magnetic field components of Jupiter
    Br, Btheta, Bphi = B_sph_components(
        r=r,
        a=RJ,
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

    # calculate Jupiter's field at callisto
    magfJ_vectors = []
    for angle in rot_theta:

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

    magfJ_vectors = np.array(magfJ_vectors) - np.mean(magfJ_vectors, axis=0)

    return magfJ_vectors


def Bext_eccentric_variation(positions):
    """
    Calculates the external field vectors with variation due to the eccentric orbit
    :param positions: array of callisto position vectors

    :return: array of external field vectors
    """

    # calculate Jupiter's field at callisto
    magfJ_vectors = []
    for pos in positions:

        # convert rotated cartesian mesh to spherical coordinates
        r, theta, phi = cartesian_to_spherical(pos)

        # Calculate legendre polynomials over a range of theta and create dictionary
        leg, dleg = legendre.quadrupole(theta)

        # calculate magnetic field components of Jupiter
        Br, Btheta, Bphi = B_sph_components(
            r=r,
            a=RJ,
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
        magfJ_vectors.append(Bcart)

    magfJ_vectors = np.array(magfJ_vectors) - np.mean(magfJ_vectors, axis=0)

    return magfJ_vectors


def Bext_full(positions, times):
    """
    Calculates the external field vectors due to both eccentric orbit and synodic rotation of Jupiter.
    :param times:
    :param positions: array of position vectors

    :return:
    """
    # time intervals and angles
    rot_theta = omega * times

    # calculate Jupiter's field at callisto
    magfJ_vectors = []
    for angle, pos in zip(rot_theta, positions):

        # convert rotated cartesian mesh to spherical coordinates
        r, theta, phi = cartesian_to_spherical(pos)

        # Calculate legendre polynomials over a range of theta and create dictionary
        leg, dleg = legendre.quadrupole(theta)

        # calculate magnetic field components of Jupiter
        Br, Btheta, Bphi = B_sph_components(
            r=r,
            a=RJ,
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