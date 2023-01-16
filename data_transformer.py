import os
import numpy as np

# # Define path to directory on any machine
# path = os.getcwd() + "\\data\\"
#
# datapoints = np.genfromtxt('data/cochlea_spiral_points.txt', dtype='float')


# Transform data points in physical space to voxels
def transform_data(data, img):
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
