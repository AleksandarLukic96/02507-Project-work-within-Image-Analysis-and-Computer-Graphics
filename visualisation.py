import math
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
from skimage import io


img = nib.load('data/cochlea.nii')
img_data = img.get_fdata()
print(img_data.shape)

img_seg = nib.load('data/cochlea_segmentation.nii')
seg_data = img_seg.get_fdata()


def show_slices(slices):
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap="gray", origin="lower")


slice_0 = img_data[449, :, :]
v_min = slice_0.min()
v_max = slice_0.max()
# slice_1 = img_data[:, 60, :]
# slice_2 = img_data[:, :, 129]
# show_slices([slice_0, slice_1, slice_2])
# plt.suptitle("Center slices for 3D image")
# plt.show()

voxel_size = 0.02449999935925007
slice_size = 200

# auxilliary function
def calculate_point_on_circle(center, radius, angle, vector1, vector2):
    result = np.array([0, 0, 0])
    for i in range(0, 3):
        result[i] = center[i] + radius * math.cos(angle) * vector1[i] + radius * math.sin(angle) * vector2[i]
    return result


def extract_points_in_image(norm_vect, data_point):
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
                visualised_slice[r, c] = img_data[rounded_point[0], rounded_point[1], rounded_point[2]]
            except IndexError:
                visualised_slice[r, c] = 0
    return visualised_slice


one_data_point = np.array([511, 216, 108])
another_data_point = np.array([514, 209, 106])
norm_vect_test = another_data_point - one_data_point


print(f"norm vector for testing: {norm_vect_test}")

testslice = extract_points_in_image(norm_vect_test, [511, 216, 108])
# useful arguments for imshow: cmap='gray', vmin=v_min, vmax=v_max
io.imshow(testslice)
io.show()



