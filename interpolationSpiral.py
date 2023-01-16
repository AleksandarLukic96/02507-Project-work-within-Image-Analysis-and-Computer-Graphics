import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import nibabel as nib

# Loading images and data
Path = 'C:\\Users\\victo\\02507_Project_ImageAnalysis_ComputerGraphics\\'
path = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea_spiral_points.txt'
pathImg = 'C:\\Users\\victo\\Dropbox\\Uni\\DTU\\MMC\\Project Image Analysis\\cochlea.nii'

img = nib.load(pathImg)

# Converting points to voxel domain
pointsMm = np.loadtxt(path)
voxelSize = img.header.get_zooms()[0] #mm
correction = np.ones(pointsMm.shape)*np.divide(img.shape,2)
pointsVoxel = np.divide(pointsMm, voxelSize) + correction

# Find the B-spline representation or its derivatives
# tck: If a tuple, a sequence of length 3, containing knots coefficients and degree of spline
# u: An array of the values of the parameter
tck, u = interpolate.splprep([pointsVoxel[:,0],pointsVoxel[:,1],pointsVoxel[:,2]], s=0)
uNew = np.linspace(0,1,10*pointsVoxel.shape[0])
xNew, yNew, zNew = interpolate.splev(uNew, tck)

# Saving interpolated points to txt file
interpolatedPoints = np.transpose(np.matrix([xNew, yNew, zNew]))
with open(Path + 'interpolatedPoints.txt','wb') as f:
    for line in interpolatedPoints:
        np.savetxt(f,line,fmt = '%.2f')

# Displaying interpolation
fig = plt.figure(2)
ax3d = fig.add_subplot(111, projection='3d')
ax3d.plot(pointsVoxel[:,0], pointsVoxel[:,1], pointsVoxel[:,2], 'r*')
ax3d.plot(xNew, yNew, zNew, 'b.', alpha = 0.5)
fig.show()
plt.show()

# Distance in between interpolated points
idx = 2
point2 = np.array((xNew[idx],yNew[idx],zNew[idx]))
point1 = np.array((xNew[idx-1],yNew[idx-1],zNew[idx-1]))
equiDistance = np.linalg.norm(point1-point2)
print(equiDistance)