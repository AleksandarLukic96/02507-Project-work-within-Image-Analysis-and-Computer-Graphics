import os
import numpy as np

width_in_mm = 14.7
height_in_mm = 9.8
depth_in_mm = 8.575

voxel_size = 0.02449999935925007

# Define path to directory on any machine
path = os.getcwd() + "\\data\\"
#print("path_to_nii: " + path_to_nii)

trn = open(path + "cochlea_spiral_points_transformed.txt", 'w')

data = np.genfromtxt('data/cochlea_spiral_points.txt', dtype='float')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

for i in range(0, x.size):
    trn.write(str(round((x[i] + width_in_mm / 2) / voxel_size)))
    trn.write(' ')
    trn.write(str(round((y[i] + height_in_mm / 2) / voxel_size)))
    trn.write(' ')
    trn.write(str(round((z[i] + depth_in_mm / 2) / voxel_size)))
    trn.write('\n')

trn.close()
