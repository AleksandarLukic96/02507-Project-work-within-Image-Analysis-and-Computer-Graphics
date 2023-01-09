import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True

x = np.linspace(-3, 3, 3)
y = np.linspace(-3, 3, 3)

x, y = np.meshgrid(x, y)

plane_equation = 0.12 * x + 0.01 * y + 1.09

fig = plt.figure()

ax = plt.axes(projection='3d')

ax.plot_surface(x, y, plane_equation, color='red')

ax.set_xlim(-10, 10) 
ax.set_ylim(-10, 10) 
ax.set_zlim(-10, 10)


plt.show()
