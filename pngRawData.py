import os
import numpy as np
from matplotlib import pyplot as plt

# Define path to directory on any machine
path = os.getcwd() + "\\data\\"
img = plt.imread(path + "saved_slices\\slice0.png")
img = plt.imread(path + "img_test.png")



num = -1
if num < 0:
    print(img.shape)
    print("img[0]:")
    print(img[0])
    print("len(img[0]):")
    print(len(img[0]))
    print("type(img): " + str(type(img)))
    print("img.shape[0]: " + str(img.shape[0]))
    print("img.shape[1]: " + str(img.shape[1]))
    

if num < -1:
    plt.imshow(img)

    x, y = np.mgrid[0:img.shape[0], 0:img.shape[1]]

    ax = plt.axes(projection="3d")
    ax.plot_surface(x, y, x*y, rstride=1, cstride=1, facecolors = img)
    plt.show()

