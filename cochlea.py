import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
from mpl_toolkits import mplot3d

path = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea_spiral_points.txt'
pathImg = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea.nii'
pathSeg = 'C:\\Users\\victo\\Dro    pbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea_segmentation.nii'

img = nib.load(pathImg)
data = img.get_fdata()

imgSegmentation = nib.load(pathSeg)
seg = imgSegmentation.get_fdata()

# Displaying slice
plt.title("Slice")
plt.xlabel("width pixel scaling")
plt.ylabel("height pixels scaling")

plt.imshow(data[:,:,2],cmap='gray')
plt.show()

# Making points from physical coordinates to voxel coordinates
pointsMm = np.loadtxt(path)

voxelSize = 0.02449999935925007
width = 600
height = 400
slices = 350

correction = np.ones(pointsMm.shape)
correction[:,0] = width/2
correction[:,1] = height/2
correction[:,2] = slices/2

pointsVoxel = np.divide(pointsMm, voxelSize) + correction


# Displaying points
# Creating figure
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection="3d")

# Creating plot
ax.scatter3D(pointsVoxel[:,0], pointsVoxel[:,1], pointsVoxel[:,2], color="green")
plt.title("Spiral points in cochlea")
plt.xlabel("width")
plt.ylabel("height")

# show plot
plt.show()

