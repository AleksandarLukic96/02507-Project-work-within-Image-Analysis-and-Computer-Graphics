import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#defining x and y arrays of the initial data set
x = np.linspace(0, 100,10)
y= 3*x**2 - np.exp(0.1*x)

# x array that will be used for interpolating new point values
x_new = np.linspace(0, 100, 100)

kind = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']

fig = plt.figure()
ax = fig.subplots()

for i in kind:
    #interpolation step
    f = interpolate.interp1d(x, y, kind = i)
    #y array that contains the interpolated data points
    y_interp = f(x_new)
    ax.plot(x_new, y_interp, alpha = 0.5, label = i)
ax.scatter(x,y)
plt.legend()
plt.show()

for i in kind:
    #interpolation step
    f = interpolate.interp1d(x, y, kind = i)
    #y array that contains the interpolated data points
    y_interp = f(x)
    ax.plot(x, y_interp, alpha = 0.5, label = i)
ax.scatter(x,y)
plt.legend()
plt.show()