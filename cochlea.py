import os
import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
from mpl_toolkits import mplot3d

# Define path to directory on any machine
path = os.getcwd() + "\\02507\\data\\cochlea\\"
#print("path_to_nii: " + path_to_nii)

path_points = path + "cochlea_spiral_points.txt"
path_normal = path + "normal_vectors_from_points.txt"
path_img = path + "cochlea.nii"
path_seg = path + "cochlea_segmentation.nii"

img = nib.load(path_img)
data = img.get_fdata()

# Get dimensions of the .nii image
epi_shape = data.shape

imgSegmentation = nib.load(path_seg)
seg = imgSegmentation.get_fdata()

# Displaying slice
plt.title("Slice")
plt.xlabel("width pixel scaling")
plt.ylabel("height pixels scaling")

plt.imshow(data[:,:,2], cmap = 'gray')
#plt.show()

# Making points from physical coordinates to voxel coordinates
points_spiral = np.loadtxt(path_points)
points_normal = np.loadtxt(path_normal)

voxelSize = 0.02449999935925007
width = epi_shape[0] #600
height = epi_shape[1] #400
slices = epi_shape[2] #350

correction_points = np.ones(points_spiral.shape)
correction_points[:,0] = width/2
correction_points[:,1] = height/2
correction_points[:,2] = slices/2

points_voxel = np.divide(points_spiral, voxelSize) + correction_points
points_normal_voxel = np.divide(points_normal, voxelSize)

# Displaying points
# Creating figure
fig = plt.figure(figsize=(10, 7))

# Set Axes accordingly to img dimensions
ax = plt.axes(projection="3d")
ax.set_xlim(0, width) 
ax.set_ylim(0, height) 
ax.set_zlim(0, slices)

# Creating plot of spiral points
ax.scatter3D(points_voxel[:,0], points_voxel[:,1], points_voxel[:,2], color="green")

# Creating plot of normal vectors:
ax.quiver(points_voxel[:,0], points_voxel[:,1], points_voxel[:,2], points_normal_voxel[:,0], points_normal_voxel[:,1], points_normal_voxel[:,2], color = "red")

# Test of normal vector in voxels
#ax.quiver(511, 216, 108, 3, -7, -2, color = "red")

#######

x = np.linspace(10, 30, 3)
y = np.linspace(10, 30, 3)

x, y = np.meshgrid(x, y)

plane_equation = 0.12 * x + 0.01 * y + 1.09

ax.plot_surface(x, y, plane_equation, color='red')

######

plt.title("Spiral points in cochlea")
plt.xlabel("width")     # x-axis = scanner-left/right 
plt.ylabel("height")    # y-axis = scanner-floor/ceiling axis
                        # z-axis = scanner-bore axis

# show plot
plt.show()

