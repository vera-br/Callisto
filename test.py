import numpy as np

R_C = 2410.3 * 1e3 # m

# define points on sphere surface
def points_on_sphere(num_points, radius=R_C):
    # Generate random spherical coordinates
    theta = np.random.uniform(0, 2 * np.pi, num_points)
    phi = np.arccos(2 * np.random.uniform(0, 1, num_points) - 1)

    # Convert spherical coordinates to Cartesian coordinates
    x = radius * np.sin(phi) * np.cos(theta)
    y = radius * np.sin(phi) * np.sin(theta)
    z = radius * np.cos(phi)

    # Return the points as a NumPy array
    return x, y, z

# create time array for one jupiter rotation (synodic)
time = np.arange(0, 10.1*3600, 60)

x, y, z = points_on_sphere(1000, R_C)

domain = np.empty((1,4))
for i in range(len(time)):
    t = np.full((len(x)), fill_value=time[i])
    domain = np.concatenate((domain, np.column_stack((t,x,y,z))))
domain = np.delete(domain, (0), axis=0) #delete first row (created with np.empty)

print(np.shape(domain))
