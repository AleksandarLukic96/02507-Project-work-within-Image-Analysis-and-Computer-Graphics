import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

pathImgInterpolated = 'data/test_nifty3.nii'
pathImg = 'data/cochlea.nii'

img = nib.load(pathImg)
imgInterpolated = nib.load(pathImgInterpolated)

img_data = img.get_fdata()
img_data_interpolated = imgInterpolated.get_fdata()

# Get voxel size from the header of the original image
voxelSize = img.header.get_zooms()[0]

# Get dimensions of the interpolated .nii image
imgShape = img_data_interpolated.shape

# Get the length of the cochlea in mm
cochleaLength = imgShape[0] * voxelSize

def getListOfDataPoints() -> list:
    """Returns a list of points in the cochlea in mm"""
    # Number of points found in the cochlea, with a distance of voxel size
    # between the points on the cochlea
    numberOfPoints = int(cochleaLength / voxelSize)
    
    # List of points in the cochlea
    points = []
    for i in range(numberOfPoints + 1):
        points.append(i * voxelSize)
    
    return points

def greenwoodFunctionForPoint(cochlealength: float, point: float) -> float:
    """Returns the frequency of a specific point in the cochlea.
    The point is in mm along the length of the cochlea"""
    x = (cochlealength - point) / cochlealength
    frequency = 165.4 * (10 ** (2.1 * x) - 0.88)
    return frequency

def greenwoodFunctionFromPoints(cochlealength: float, points: list) -> list:
    """Returns a list of frequencies that can be heard 
    at the given points in the cochlea according to its length"""
    
    frequencies = []
    for point in points:
        x = (cochlealength - point) / cochlealength
        frequency = 165.4 * (10 ** (2.1 * x) - 0.88)
        frequencies.append(frequency)
    
    return frequencies

# Generate color list from frequencies
def getFrequencyColors(frequencies: list):
    colors = []
    for frequency in frequencies:
        if frequency < 125.0:
            colors.append('red')
        elif frequency < 250.0:
            colors.append('orange')
        elif frequency < 500.0:
            colors.append('yellow')
        elif frequency < 1000.0:
            colors.append('green')
        elif frequency < 2000.0:
            colors.append('lime')
        elif frequency < 4000.0:
            colors.append('lightblue')
        elif frequency < 8000.0:
            colors.append('blue')
        elif frequency < 16000.0:
            colors.append('violet')
        else:
            colors.append('purple')
    return colors

# Get the list of data points in mm from the base of the cochlea
dataPoints = getListOfDataPoints()

# Get list of frequencies
frequencyList = greenwoodFunctionFromPoints(cochleaLength, dataPoints)

# Get colors for frequencies
colors = getFrequencyColors(frequencyList)

# Data points to plot
data = np.genfromtxt('interpolatedPoints.txt', dtype = 'float')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Displaying points
# Creating figure
fig = plt.figure(figsize=(10, 7))

# Set Axes accordingly to img dimensions
ax = plt.axes(projection="3d")

# Creating plot of spiral points
ax.scatter3D(x, y, z, c=colors)

plt.title("Spiral points in cochlea")
plt.xlabel("width")     # x-axis = scanner-left/right 
plt.ylabel("height")    # y-axis = scanner-floor/ceiling axis
                        # z-axis = scanner-bore axis
                        
# show plot
plt.show()