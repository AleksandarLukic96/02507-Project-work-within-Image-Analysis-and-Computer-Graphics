import os
import sys
import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from nibabel.affines import apply_affine
from skimage import io, segmentation, measure
from skimage.morphology import dilation, erosion, disk

slice_size = 200


# auxilliary function - deprecated
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


