import numpy as np

# Function to calculate the frequency that can be heard 
# at certain points in the cochlea according to its length
def greenwoodFunction(cochlealength):
    # Number of points to find the frequency for, with a distance of
    # 25 microns between the points on the cochlea
    numPoints = round(cochlealength * 1000 / 25)
    
    frequencies = []
    for point in np.linspace(0, cochlealength, numPoints, endpoint=True):
        x = (cochlealength - point) / cochlealength
        frequency = 165.4 * (10 ** (2.1 * x) - 0.88)
        frequencies.append(frequency)
    
    return frequencies