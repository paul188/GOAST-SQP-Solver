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

positions = []
normals = []

# Open the file and read it line by line
with open('normals.txt', 'r') as file:
    n = sum(1 for _ in file)

    # Reset file pointer to the beginning
    file.seek(0)

    # Read the first n lines as positions
    for _ in range(n // 2):  # You should know n, the number of position vectors
        position_line = file.readline().strip()
        #numbers = [float(x) for x in position_line.split()]
        positions.append([float(x) for x in position_line.split()])

    # Read the second n lines as normal vectors
    for _ in range(n // 2):  # The same n for the normals
        normal_line = file.readline().strip()
        normals.append([float(x) for x in normal_line.split()])

# Extract x, y, z components separately
x_coords = [normal[0] for normal in normals]
y_coords = [normal[1] for normal in normals]
z_coords = [normal[2] for normal in normals]

print("X test: ")
print(x_coords)

print("Y test: ")
print(y_coords)

print("Z test: ")
print(z_coords)

ax.scatter(x_coords, y_coords, z_coords, color='b',s=0.001)  # Scatter plot of points on the sphere

# Labels and aspect ratio
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_box_aspect([1, 1, 1])  # Ensure proper scaling

ax.set_title('Plot of the Gauss image of the surface on the unit sphere')

plt.show()
