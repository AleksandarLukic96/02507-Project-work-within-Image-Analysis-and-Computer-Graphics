import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt
import os

# Define path to directory on any machine
path = os.getcwd() + "\\02507\\data\\cochlea\\"
#print("path: " + path)

# Choose .nii file for examination 
nib_file = "cochlea.nii"

# Load .nii image and its data
epi_img = nib.load(path + nib_file)
epi_img_data = epi_img.get_fdata()

# Prints the dimensions of the .nii image
epi_shape = epi_img_data.shape
#print(epi_shape)

# Voxel size definition
voxel_size = 0.02449999935925007

# Values collected from via shape of nii image
width_in_mm = epi_shape[0] * voxel_size                      
height_in_mm = epi_shape[1] * voxel_size
depth_in_mm = epi_shape[2] * voxel_size
#print("(" + str(width_in_mm) + ", " + str(height_in_mm) + ", " + str(depth_in_mm) + ")")

# Find vector between two points
def vector3d_between_two_points(m1, m2):
    return np.array([m2[0] - m1[0], m2[1] - m1[1], m2[2] - m1[2]])

file = open(path + 'normal_vectors_from_points.txt', 'w')

data = np.genfromtxt(path + 'cochlea_spiral_points.txt', dtype = 'float')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

for i in range(0, x.size-1):
    v1 = np.array([x[i], y[i], z[i]])
    v2 = np.array([x[i+1], y[i+1], z[i+1]])
    v3 = vector3d_between_two_points(v1, v2)
    file.write(str(v3[0]))
    file.write(' ')
    file.write(str(v3[1]))
    file.write(' ')
    file.write(str(v3[2]))
    file.write('\n')
    
file.close()