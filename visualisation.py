import os
import sys
import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage import io, measure
from skimage.morphology import dilation, erosion, disk
from interpolationVoxelValues import trilinearInterpolation

# ########## GLOBAL VARIABLES #########
# Currently the size of the slices to draw is hard-coded to be this number:
slice_size = 200

# This indicates how thick a border we want when removing redundant data
border_thickness = 2
#######################################

# Define path to directory on any machine
path = os.getcwd() + "\\data\\"

# Loading some of the relevant data
img = nib.load(path + "cochlea.nii")
img_data = img.get_fdata()

img_seg = nib.load(path + "cochlea_segmentation.nii")
seg_data = img_seg.get_fdata()

new_data_points = np.genfromtxt("interpolatedPoints.txt")
new_vectors = np.genfromtxt("normalVec.txt")


# This function uses a scatterplot to plot the shape of the cochlea in an
# interactive 3D format. Give it a dataset as input to also plot that
def plot_3d_of_cochlea(resolution=10, data=None):
    xyzdata = []
    for x in range(0, 600, resolution):
        for y in range(0, 400, resolution):
            for z in range(0, 350, resolution):
                if seg_data[x, y, z] == 1:
                    xyzdata.append([x, y, z])
    xyzdata = np.array(xyzdata)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    zdata = xyzdata[:, 2]
    ax.scatter3D(xyzdata[:, 0], xyzdata[:, 1], zdata, c=zdata, alpha=0.1)
    if data is not None:
        ax.scatter3D(data[:, 0], data[:, 1], data[:, 2], cmap='Reds')
    plt.show()


# Extract data from the 3D image and draw it as a 2D slice
def extract_points_in_image(norm_vect, data_point, input_img, segmentation=False, mask=None):
    # This method uses the much simpler approach suggested by Martin, in hopes of achieving
    # a faster execution time
    visualised_slice = np.zeros((slice_size, slice_size))
    i_vect = np.cross(norm_vect, np.array([0, 0, 1]))
    if all(v == 0 for v in i_vect):
        i_vect = norm_vect
    i_vect_n = i_vect / np.linalg.norm(i_vect)
    j_vect = np.cross(norm_vect, i_vect)
    j_vect_n = j_vect / np.linalg.norm(j_vect)

    for r in range(0, slice_size):
        for c in range(0, slice_size):
            if mask is not None and mask[r, c] == 0:
                continue
            point_in_data = (r - slice_size / 2) * i_vect_n + (c - slice_size / 2) * j_vect_n + data_point
            try:
                if segmentation:
                    visualised_slice[r, c] = input_img[int(round(point_in_data[0])), int(round(point_in_data[1])), int(round(point_in_data[2]))]
                else:
                    visualised_slice[r, c] = trilinearInterpolation(point_in_data, input_img)
            except IndexError:
                visualised_slice[r, c] = 0
    return visualised_slice


# This function extracts a slice from the given plane in a segmentation image
# and uses BLOB analysis techniques to create a mask
def generate_mask(norm_vect, data_point, input_img):
    seg_slice = extract_points_in_image(norm_vect, data_point, input_img, segmentation=True)
    seg_slice = erosion(seg_slice, disk(4))

    label_slice = measure.label(seg_slice, connectivity=1)
    region_props = measure.regionprops(label_slice)

    label_slice_filter = label_slice

    center = np.array([int(slice_size / 2), int(slice_size / 2)])
    min_dist = slice_size
    min_dist_label = 0
    for region in region_props:
        distance = math.dist(center, region.centroid)
        if distance < min_dist:
            min_dist = distance
            min_dist_label = region.label

    label_slice_filter[label_slice_filter != min_dist_label] = 0

    label_slice_filter = dilation(label_slice_filter, disk(12 + border_thickness))
    label_slice_filter = erosion(label_slice_filter, disk(10))

    return label_slice_filter


# Using the functions above, draw slices for all available data points and save
# it as a combined nifty file
def draw_nifty(filename, affine):
    list_of_slices = []
    no_slices = new_data_points.shape[0]
    for i in range(0, no_slices - 1): # 0, no_slices - 1
        point = new_data_points[i]
        vect = new_vectors[i]

        mask = generate_mask(vect, point, seg_data)
        test_slice = extract_points_in_image(vect, point, img_data, mask=mask)

        cutout = test_slice
        cutout[mask == 0] = 0

        list_of_slices.append(cutout)
        # write progress in the terminal
        sys.stdout.write("\rCompleted slice %i out of %i" % (i+1, no_slices))
        sys.stdout.flush()

    array_of_slices = np.array(list_of_slices)
    nifti_file = nib.Nifti1Image(array_of_slices, affine)
    nib.save(nifti_file, path + filename)


test_affine = np.array([[0.0245, 0, 0, 0],
                        [0, 0.0245, 0, 0],
                        [0, 0, 0.0245, 0],
                        [0, 0, 0, 1]])


draw_nifty("test_nifty3.nii", test_affine)

