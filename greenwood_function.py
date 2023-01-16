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