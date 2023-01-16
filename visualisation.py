import os
import sys
import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from nibabel.affines import apply_affine
from skimage import io, segmentation, measure
from skimage.morphology import dilation, erosion, disk
from skimage.color import label2rgb


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
print(img_data.shape)

img_seg = nib.load(path + "cochlea_segmentation.nii")
seg_data = img_seg.get_fdata()

data_points = np.genfromtxt(path + "cochlea_spiral_points.txt", dtype='float')


# Transform data points in physical space to voxels
def transform_data(data):
    trn = open('data/cochlea_spiral_points_transformed.txt', 'w')

    voxel_size = img.header['pixdim'][1]
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    for i in range(0, x.size):
        trn.write(str(round((x[i] + img.header['qoffset_x']) / voxel_size)))
        trn.write(' ')
        trn.write(str(round((y[i] + img.header['qoffset_y']) / voxel_size)))
        trn.write(' ')
        trn.write(str(round((z[i] - img.header['qoffset_z']) / voxel_size)))
        trn.write('\n')

    trn.close()


# Load the transformed data after transformation
transform_data(data_points)
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


# auxilliary function
def calculate_point_on_circle(center, radius, angle, vector1, vector2):
    result = np.array([0, 0, 0])
    for i in range(0, 3):
        result[i] = center[i] + radius * math.cos(angle) * vector1[i] + radius * math.sin(angle) * vector2[i]
    return result


# This function draws the slices based on a normal vector and a central data point.
# It also needs to know which image to extract data from
def extract_points_in_image_old(norm_vect, data_point, input_img):
    # initialise the image in which we'll place our datapoints:
    visualised_slice = np.zeros((slice_size, slice_size))
    # find a vector on the plane parallel to the x-axis:
    vert_vect = np.cross(norm_vect, np.array([0, 1, 0]))
    vert_vect_n = vert_vect/np.linalg.norm(vert_vect)
    # find another vector on the plane, that is not parallel to vert_vect
    w_vect = np.cross(norm_vect, vert_vect)
    w_vect_n = w_vect/np.linalg.norm(w_vect)

    for r in range(0, slice_size):
        for c in range(0, slice_size):
            # define a vector between the data point and the current pixel to draw
            x = c - slice_size / 2
            y = r - slice_size / 2
            point_vect = np.array([x, y])
            # calculate the length of this vector
            radius = np.linalg.norm(point_vect)
            # calculate the angle between the previous vector and the vertical vector
            if radius == 0:
                angle = 0
            else:
                simple_vert = np.array([1, 0])
                angle = math.acos(simple_vert.dot(point_vect) / (np.linalg.norm(simple_vert) * radius))
            # the angle found is always the smallest angle between the vectors. So, if we've reached
            # the bottom part of the image, it will find the wrong angle. Thus, the conversion below is needed.
            if point_vect[1] > 0:
                angle = math.pi * 2 - angle

            point_in_data = calculate_point_on_circle(data_point, radius, angle, vert_vect_n, w_vect_n)
            # now a value is fetched simply from the nearest point in the 3D image. When interpolation
            # is implemented, it can be used here
            rounded_point = np.array([int(i) for i in point_in_data])
            try:
                visualised_slice[r, c] = input_img[rounded_point[0], rounded_point[1], rounded_point[2]]
            except IndexError:
                visualised_slice[r, c] = 0
    return visualised_slice


def extract_points_in_image(norm_vect, data_point, input_img):
    # This method uses the much simpler approach suggested by Martin, in hopes of achieving
    # a faster execution time
    visualised_slice = np.zeros((slice_size, slice_size))
    i_vect = np.cross(norm_vect, np.array([0, 0, 1]))
#   print(f"i_vect before: {i_vect}")
    if all(v == 0 for v in i_vect):
        i_vect = norm_vect
#    print(f"i_vect after: {i_vect}")
    i_vect_n = i_vect/np.linalg.norm(i_vect)
    j_vect = np.cross(norm_vect, i_vect)
#    print(f"j_vect: {j_vect}")
    j_vect_n = j_vect / np.linalg.norm(j_vect)

    for r in range(0, slice_size):
        for c in range(0, slice_size):
            point_in_data = (r - slice_size/2) * i_vect_n + (c - slice_size/2) * j_vect_n + data_point
            rounded_point = np.array([int(i) for i in point_in_data])
            try:
                visualised_slice[r, c] = input_img[rounded_point[0], rounded_point[1], rounded_point[2]]
            except IndexError:
                visualised_slice[r, c] = 0
    return visualised_slice


# A bunch of rotations. Currently not used, but might come in handy later
def rotation1(angle):
    angle = angle * math.pi / 180
    return np.array([[1, 0, 0],
                     [0, math.cos(angle), -math.sin(angle)],
                     [0, math.sin(angle), math.cos(angle)]])


def rotation2(angle):
    angle = angle * math.pi / 180
    return  np.array([[math.cos(angle), 0, math.sin(angle)],
                      [0, 1, 0],
                      [-math.sin(angle), 0, math.cos(angle)]])


def rotation3(angle):
    angle = angle * math.pi / 180
    return np.array([[math.cos(angle), -math.sin(angle), 0],
                     [math.sin(angle), math.cos(angle), 0],
                     [0, 0, 1]])


def make_vectors_naive(spiral_points):
    result = np.zeros(spiral_points.shape)
    for i in range(0, spiral_points.shape[0] - 1):
        result[i] = spiral_points[i+1] - spiral_points[i]
    return result


def draw_nifty(filename, affine):
    list_of_slices = []
    no_slices = transformed_data.shape[0]
    for i in range(0, no_slices-1):
        one_point = transformed_data[i]
        another_point = transformed_data[i+1]
        test_vect = another_point - one_point

        test_slice = extract_points_in_image(test_vect, one_point, img_data)
        test_slice_seg = extract_points_in_image(test_vect, one_point, seg_data)
        test_slice_seg = erosion(test_slice_seg, disk(2))

        label_slice = measure.label(test_slice_seg, connectivity=1)
        region_props = measure.regionprops(label_slice)

        label_slice_filter = label_slice

        center = np.array([int(slice_size/2), int(slice_size/2)])
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


draw_nifty("test_nifty2.nii", np.eye(4))





