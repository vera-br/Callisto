import numpy as np
from scipy.spatial.transform import Rotation as ROT

"""
For maths formulas and coordinate transformations.
"""


def threeDmag(A):
    return np.sqrt(A[0] ** 2 + A[1] ** 2 + A[2] ** 2)


def unit_spherical_to_cartesian(a, theta, phi):
    """
    Convert from 3D spherical coordinates to cartesian coordinates
    Jacobian found: https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
    :param a: 3D np.ndarray of the form [A_radial, A_polar, A_azimuthal]
    :param theta: polar angle
    :param phi: azimuthal angle
    :return:
    """

    # x component
    Bx = (
        np.sin(theta) * np.cos(phi) * a[0]
        + np.cos(theta) * np.cos(phi) * a[1]
        - np.sin(phi) * a[2]
    )

    # y component
    By = (
        np.sin(theta) * np.sin(phi) * a[0]
        + np.cos(theta) * np.sin(phi) * a[1]
        + np.cos(phi) * a[2]
    )

    # z component
    Bz = np.cos(theta) * a[0] - np.sin(theta) * a[1]

    return np.array([Bx, By, Bz])


def unit_cartesian_to_spherical(a, theta, phi):
    """
    Convert from 3D cartesian coordinates to spherical coordinates
    Jacobian found: https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
    :param a: 3D np.ndarray of the form [A_radial, A_polar, A_azimuthal]
    :param theta: polar angle
    :param phi: azimuthal angle
    :return:
    """

    ax = a[0]
    ay = a[1]
    az = a[2]

    # radial component
    Br = (
        np.sin(theta) * np.cos(phi) * ax
        + np.sin(theta) * np.sin(phi) * ay
        + np.cos(theta) * az
    )

    # theta component
    Btheta = (
        np.cos(theta) * np.cos(phi) * ax
        + np.cos(theta) * np.sin(phi) * ay
        - np.sin(theta) * az
    )

    # phi component
    Bphi = -np.sin(phi) * ax + np.cos(phi) * ay

    return np.array([Br, Btheta, Bphi])


def threeDx_rotation(A, theta):
    """
    Rotate 3D vector about x axis
    :param A: 3D vector
    :param theta: rotation angle
    :return:
    rotated vector
    """
    r1 = [1, 0, 0]
    r2 = [0, np.cos(theta), -np.sin(theta)]
    r3 = [0, np.sin(theta), np.cos(theta)]

    t = np.array([r1, r2, r3])

    return np.matmul(t, A)


def threeDz_rotation(A, theta):
    """
    Rotate 3D vector about z axis
    :param A: 3D vector
    :param theta: rotation angle
    :return:
    rotated vector
    """

    zero = np.full_like(theta, 0)
    one = np.full_like(theta, 1)
    r1 = [np.cos(theta), -np.sin(theta), 0]
    r2 = [np.sin(theta), np.cos(theta), 0]
    r3 = [0, 0, 1]

    t = np.array([r1, r2, r3])

    return np.matmul(t, A)


def spherical_to_cartesian(A):
    """
    Converts spherical vector to cartesian vector
    :param A:
    :return:
    3D vector
    """
    x = A[0] * np.sin(A[1]) * np.cos(A[2])
    y = A[0] * np.sin(A[1]) * np.sin(A[2])
    z = A[0] * np.cos(A[1])

    return np.array([x, y, z])


def cartesian_to_spherical(A):
    """
    Converts cartesian vector to spherical vector
    :param A: cartesian vector
    :return: spherical vector
    """
    r = threeDmag(A)
    theta = np.arccos(A[2] / r)  # arccos(z/r)
    phi = np.arctan2(A[1], A[0])  # arctan(y/x)

    return r, theta, phi


def cylindrical_to_cartesian(A, theta):

    r = np.linalg.norm(A)
    Ax = r * np.cos(theta)
    Ay = r * np.sin(theta)
    Az = A[-1]

    return np.array([Ax, Ay, Az])


def unit_vector(vector):
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def general_rotation(vector, rotation_angle, axis_of_rotation):

    x, y, z = vector.transpose()

    rotation = ROT.from_rotvec(rotation_angle * axis_of_rotation)

    XYZ = np.array([x.ravel(), y.ravel(), z.ravel()]).transpose()
    XYZrot = rotation.apply(XYZ)

    xrot = XYZrot[:, 0].reshape(x.shape)
    yrot = XYZrot[:, 1].reshape(y.shape)
    zrot = XYZrot[:, 2].reshape(z.shape)

    return np.array([xrot, yrot, zrot]).transpose()
