import numpy as np
np.set_printoptions(precision=7, suppress=True)
import nibabel as nib
from matplotlib import pyplot as plt

Path = 'C:\\Users\\victo\\02507_Project_ImageAnalysis_ComputerGraphics\\'
path = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea_spiral_points.txt'
pathImg = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea.nii'
pathSeg = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea_segmentation.nii'

img = nib.load(pathImg)
data = img.get_fdata()

imgSegmentation = nib.load(pathSeg)
seg = imgSegmentation.get_fdata()

# Making points from physical coordinates to voxel coordinates
# pointsMm = np.loadtxt(path)
# voxelSize = img.header.get_zooms()[0] #mm
# correction = np.ones(pointsMm.shape)*np.divide(img.shape,2)
# pointsVoxel = np.divide(pointsMm, voxelSize) + correction
# roundPointsVoxel = np.around(np.divide(pointsMm, voxelSize) + correction)

# Interpolated points in voxels
pointsVoxel = np.loadtxt(Path + 'interpolatedPoints.txt')

# Finding circle center and radius from 3 points
# https://itecnote.com/tecnote/python-find-arc-circle-equation-given-three-points-in-space-3d/
normalVec = np.zeros(shape=(pointsVoxel.shape), dtype='object')

for n, point in enumerate(pointsVoxel):
    if n == (pointsVoxel.shape[0]-2):
        break

    idx = n + 1
    p1 = pointsVoxel[n, :]
    p2 = pointsVoxel[idx, :]
    p3 = pointsVoxel[idx+1, :]

    p1p2 = np.subtract(p2, p1)
    p1p3 = np.subtract(p3, p1)
    p2p3 = np.subtract(p3, p2)

    if np.all(np.cross(p1p2, p1p3) <= 3.0e-12):
        normalVec[idx, :] = p1p2
        continue

    a = np.linalg.norm(p2p3)
    b = np.linalg.norm(p1p3)
    c = np.linalg.norm(p1p2)

    s = (a + b + c) / 2
    R = a * b * c / 4 / np.sqrt(s * (s - a) * (s - b) * (s - c))
    b1 = a * a * (b * b + c * c - a * a)
    b2 = b * b * (a * a + c * c - b * b)
    b3 = c * c * (a * a + b * b - c * c)
    center = np.column_stack((p1, p2, p3)).dot(np.hstack((b1, b2, b3)))
    center /= b1 + b2 + b3

    p2c = np.subtract(center, p2)
    n = np.cross(p1p2, p1p3)
    nVec = np.cross(p2c, n)
    normalVec[idx,:] = nVec / np.linalg.norm(nVec) * 20

normalVec[0,:] = np.subtract(pointsVoxel[1,:],pointsVoxel[0,:])# first normal vector
normalVec[normalVec.shape[0]-1,:] = normalVec[(normalVec.shape[0]-2),:]# last normal vector

# Saving txt file with normal vectors
normalVecMatrix = np.matrix(normalVec)
with open(Path + 'normalVec.txt', 'wb') as f:
    for line in normalVecMatrix:
        np.savetxt(f, line, fmt='%.2f')

normalVecPlot = np.hstack((pointsVoxel,normalVec))

# Displaying points
# Creating figure
fig = plt.figure(figsize=(10, 7))
ax = plt.axes(projection="3d")

# Creating plot
ax.scatter3D(pointsVoxel[:,0], pointsVoxel[:,1], pointsVoxel[:,2], color="green")
ax.quiver(normalVecPlot[1:(normalVecPlot.shape[0]-1),0],normalVecPlot[1:(normalVecPlot.shape[0]-1),1],
          normalVecPlot[1:(normalVecPlot.shape[0]-1),2],normalVecPlot[1:(normalVecPlot.shape[0]-1),3],
          normalVecPlot[1:(normalVecPlot.shape[0]-1),4],normalVecPlot[1:(normalVecPlot.shape[0]-1),5])
plt.title("Spiral points in cochlea")
plt.xlabel("width")
plt.ylabel("height")

# show plot
plt.show()