import os
import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
from mpl_toolkits import mplot3d
from skimage import img_as_ubyte
from skimage.color import rgb2gray
from PIL import Image
import matplotlib.image as mpimg


# Define path to directory on any machine
path = os.getcwd() + "\\data\\"
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

################################################################################
img = Image.open(path + "img_test.png")
a = img.convert("P", palette=Image.ADAPTIVE, colors=8)
a.save(path + "img_test2.png")

img2 = plt.imread(path + "img_test2.png")

print("type(img2): " + str(type(img2)))
print("img2.shape[0]: " + str(img2.shape[0]))
print("img2.shape[1]: " + str(img2.shape[1]))
print("img2.shape[2]: " + str(img2.shape[2]))

print(img2)
plt.imshow(img2)

x, y = np.mgrid[0:img2.shape[0], 0:img2.shape[1]]

ax = plt.axes(projection="3d")
ax.plot_surface(x, y, x*y, rstride=10, cstride=10, facecolors = img2)
plt.show()
################################################################################

################################################################################
img = Image.open(path + "test_slice.png")
a = img.convert("P", palette=Image.ADAPTIVE, colors=8)
a.save(path + "test_slice2.png")

#img1 = img_as_ubyte(plt.imread(path + "test_slice2.png"))
img1 = plt.imread(path + "test_slice2.png")

print("type(img1): " + str(type(img1)))
print("img1.shape[0]: " + str(img1.shape[0]))
print("img1.shape[1]: " + str(img1.shape[1]))

print(img1)
plt.imshow(img1)

x, y = np.mgrid[0:img1.shape[0], 0:img1.shape[1]]

ax = plt.axes(projection="3d")
ax.plot_surface(x, y, x*y, rstride=1, cstride=1, facecolors = img1)
plt.show()
################################################################################

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
ax.quiver(points_voxel[:,0], points_voxel[:,1], points_voxel[:,2], points_normal_voxel[:,0], points_normal_voxel[:,1], points_normal_voxel[:,2], color = "blue")

# Test of normal vector in voxels
#ax.quiver(511, 216, 108, 3, -7, -2, color = "red")

#######
for i in range(len(points_voxel)-1, 0, -1):
    x0 = round(points_voxel[i,0])
    y0 = round(points_voxel[i,1])
    z0 = round(points_voxel[i,2])
    #print("x0: " + str(x0) + " y0: " + str(y0) + " z0: " + str(z0))

    x = np.linspace(x0-5, x0+5, 3)
    y = np.linspace(y0-5, y0+5, 3)
    x, y = np.meshgrid(x, y)

    a = round(points_normal_voxel[i][0])
    b = round(points_normal_voxel[i][1])
    c = round(points_normal_voxel[i][2])
    #print("a: " + str(a) + " b: " + str(b) + " c: " + str(c))

    # Rewritten equation: a*(x-x0) / -c + b*(y-y0) / -c + z0 = z
    plane_equation = (a*(x-x0) / (-c)) + (b*(y-y0) / (-c)) + z0
    #print(f"({a}*(x-{x0}) / (-{c})) + ({b}*(y-{y0}) / (-{c})) + {z0}")

    ax.plot_surface(x, y, plane_equation, color = 'red')
######

plt.title("Spiral points in cochlea")
plt.xlabel("width")     # x-axis = scanner-left/right 
plt.ylabel("height")    # y-axis = scanner-floor/ceiling axis
                        # z-axis = scanner-bore axis

# show plot
plt.show()

