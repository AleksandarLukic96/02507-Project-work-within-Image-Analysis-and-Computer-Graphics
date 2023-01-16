import os
import sys
import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage import io, segmentation, measure
from skimage.morphology import dilation, erosion, disk
from data_transformer import transform_data

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

data_points = np.genfromtxt(path + "cochlea_spiral_points.txt", dtype='float')

# Load the transformed data after transformation
transform_data(data_points, img)
transformed_data = np.genfromtxt(path + "cochlea_spiral_points_transformed.txt", dtype='float')


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
def extract_points_in_image(norm_vect, data_point, input_img):
    # This method uses the much simpler approach suggested by Martin, in hopes of achieving
    # a faster execution time
    visualised_slice = np.zeros((slice_size, slice_size))
    i_vect = np.cross(norm_vect, np.array([0, 0, 1]))
    #   print(f"i_vect before: {i_vect}")
    if all(v == 0 for v in i_vect):
        i_vect = norm_vect
    #    print(f"i_vect after: {i_vect}")
    i_vect_n = i_vect / np.linalg.norm(i_vect)
    j_vect = np.cross(norm_vect, i_vect)
    #    print(f"j_vect: {j_vect}")
    j_vect_n = j_vect / np.linalg.norm(j_vect)

    for r in range(0, slice_size):
        for c in range(0, slice_size):
            point_in_data = (r - slice_size / 2) * i_vect_n + (c - slice_size / 2) * j_vect_n + data_point
            rounded_point = np.array([int(i) for i in point_in_data])
            try:
                visualised_slice[r, c] = input_img[rounded_point[0], rounded_point[1], rounded_point[2]]
            except IndexError:
                visualised_slice[r, c] = 0
    return visualised_slice


# Using the function above, draw slices for all available data points and save
# it as a combined nifty file
def draw_nifty(filename, affine):
    list_of_slices = []
    no_slices = transformed_data.shape[0]
    for i in range(0, no_slices - 1):
        one_point = transformed_data[i]
        another_point = transformed_data[i + 1]
        test_vect = another_point - one_point

        test_slice = extract_points_in_image(test_vect, one_point, img_data)
        test_slice_seg = extract_points_in_image(test_vect, one_point, seg_data)
        test_slice_seg = erosion(test_slice_seg, disk(2))

        label_slice = measure.label(test_slice_seg, connectivity=1)
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

        label_slice_filter = dilation(label_slice_filter, disk(10 + border_thickness))
        label_slice_filter = erosion(label_slice_filter, disk(10))

        cutout = test_slice
        cutout[label_slice_filter == 0] = 0

        list_of_slices.append(cutout)
        # write progress in the terminal
        sys.stdout.write("\rCompleted slice %i out of %i" % (i, no_slices))
        sys.stdout.flush()

    array_of_slices = np.array(list_of_slices)
    nifti_file = nib.Nifti1Image(array_of_slices, affine)
    nib.save(nifti_file, path + filename)


test_affine = np.array([[0.0245, 0, 0, 0],
                        [0, 0.0245, 0, 0],
                        [0, 0, 0.0245, 0],
                        [0, 0, 0, 1]])


# draw_nifty("test_nifty2.nii", test_affine)


# Calculate the length between all the given data points (source is currently hard-coded)
def estimate_length_raw_datapoints():
    sum = 0
    for i in range(0, data_points.shape[0] - 1):
        sum = sum + math.dist(data_points[i], data_points[i + 1])
    return sum


print(f"App. expected length of cochlea when fully unfolded: {estimate_length_raw_datapoints()}")
