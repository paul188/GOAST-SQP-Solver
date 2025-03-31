from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Draw sphere
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
ax.plot_surface(x, y, z, color="g", alpha=0.3)  # Spherical surface

# Generate random points on a unit sphere
'''
theta = np.random.uniform(0, 2*np.pi, 100)  # Azimuthal angle
phi = np.arccos(2 * np.random.uniform(0, 1, 100) - 1)  # Elevation angle (uniform sampling)

# Convert to Cartesian coordinates
x_pts = np.sin(phi) * np.cos(theta)
y_pts = np.sin(phi) * np.sin(theta)
z_pts = np.cos(phi)'

# Scatter plot of points on the sphere
ax.scatter(x_pts, y_pts, z_pts, color='b')'
'''

x_coords = []
y_coords = []
z_coords = []

# Open the file and read each line
with open("/home/paul_johannssen/Desktop/masterarbeit/goast_old/examples/QuadMeshExamples/plotting/gauss_image_data_final.txt", "r") as file:
    for line in file:
        # Split the line by whitespace and extract x, y, z
        coords = line.split()
        x, y, z = float(coords[0]), float(coords[1]), float(coords[2])
        
        # Append the coordinates to respective arrays
        x_coords.append(x)
        y_coords.append(y)
        z_coords.append(z)

ax.scatter(x_coords, y_coords, z_coords, color='b')  # Scatter plot of points on the sphere

# Labels and aspect ratio
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_box_aspect([1, 1, 1])  # Ensure proper scaling

ax.set_title('Random Points on a Unit Sphere')

plt.show()
