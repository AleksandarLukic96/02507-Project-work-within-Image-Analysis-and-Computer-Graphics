import os
import sys
import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage import io, measure, img_as_ubyte
from skimage.color import gray2rgb, rgb2hsv, hsv2rgb

# Setting the path and loading the 3D file
path = os.getcwd() + "\\data\\"
nifty_data = nib.load(path + "test_nifty3.nii").get_fdata()

# Extract the prettiest cross section from the nifti
cross_sec = nifty_data[:, 100, :]

# bringing the image values down to [0, 1]
cross_sec = cross_sec - cross_sec[800, 165]
max_val = cross_sec.max()
cross_sec = cross_sec / max_val

# converting the image to HSV for easier colour processing
cross_sec_rgb = gray2rgb(cross_sec)
cross_sec_hsv = rgb2hsv(cross_sec_rgb)


# auxilliary function for showing all hsv channels separately
def show_comparison(hsv_img):
    f, ax = plt.subplots(1, 3)
    ax[0].imshow(hsv_img[:, :, 0], cmap='hsv')
    ax[0].set_title("Hue image")
    ax[1].imshow(hsv_img[:, :, 1])
    ax[1].set_title("Saturation image")
    ax[2].imshow(hsv_img[:, :, 2])
    ax[2].set_title("Value image")
    io.show()


# setting the Saturation channel to equal the Value channel (which currently holds the entire gray-scale image)
cross_sec_edit = cross_sec_hsv
cross_sec_edit[:, :, 1] = cross_sec_edit[:, :, 2] + 0.2


# Looping through the image to set the desired Hue values
i = 0
for r in range(0, cross_sec.shape[0]):
    i = i + 1/1110
    for c in range(0, cross_sec.shape[1]):
        cross_sec_edit[r, c, 0] = i


# showing the result
show_comparison(cross_sec_edit)
cross_sec_edit = hsv2rgb(cross_sec_edit)
io.imshow(cross_sec_edit)
io.show()





