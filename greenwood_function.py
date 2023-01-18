import numpy as np

def greenwoodFunction(cochlealength: float) -> list:
    """Returns a list of frequencies that can be heard 
    at certain points in the cochlea according to its length"""
    
    # Number of points to find the frequency for, with a distance of
    # 25 microns between the points on the cochlea
    numPoints = round(cochlealength * 1000 / 25)
    
    frequencies = []
    for point in np.linspace(0, cochlealength, numPoints, endpoint=True):
        x = (cochlealength - point) / cochlealength
        frequency = 165.4 * (10 ** (2.1 * x) - 0.88)
        frequencies.append(frequency)
    
    return frequencies

def greenwoodFunctionPoint(cochlealength: float, point: float) -> float:
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
    for i in range(len(frequencies)):
        if frequencies[i] < 125.0:
            colors.append('red')
        elif frequencies[i] < 250.0:
            colors.append('orange')
        elif frequencies[i] < 500.0:
            colors.append('yellow')
        elif frequencies[i] < 1000.0:
            colors.append('green')
        elif frequencies[i] < 2000.0:
            colors.append('lime')
        elif frequencies[i] < 4000.0:
            colors.append('lightblue')
        elif frequencies[i] < 8000.0:
            colors.append('blue')
        elif frequencies[i] < 16000.0:
            colors.append('violet')
        else:
            colors.append('purple')
    return colors