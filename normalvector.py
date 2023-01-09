# This file takes a .txt file of points and calculates the cooresponding normal vector for each point
# in order to define each point's plane to slice the .nii image

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

# Get dimensions of the .nii image
epi_shape = epi_img_data.shape
#print(epi_shape)

# Voxel size definition
voxel_size = 0.02449999935925007

# Values collected from via shape of nii image
width_in_mm = epi_shape[0] * voxel_size                      
height_in_mm = epi_shape[1] * voxel_size
depth_in_mm = epi_shape[2] * voxel_size

# Find vector between two points
def vector3d_between_two_points(v1, v2):
    return v2 - v1

file = open(path + 'normal_vectors_from_points.txt', 'w')

data = np.genfromtxt(path + 'cochlea_spiral_points.txt', dtype = 'float')
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Temp array to store normal vectors in reverse order
M = []

for i in range(x.size-1, 0, -1):
    #print(i)
    v1 = np.array([x[i], y[i], z[i]])
    #print(v1)
    v2 = np.array([x[i-1], y[i-1], z[i-1]])
    #print(v2)
    v3 = vector3d_between_two_points(v1, v2) * (-1)
    #print(v3)
    M.append([v3[0], v3[1], v3[2]])

print(x.size)
print(len(M))
print(str(M[x.size-2]))

# Edge case:
file.write(str(M[len(M)-1][0]))
file.write(' ')
file.write(str(M[len(M)-1][1]))
file.write(' ')
file.write(str(M[len(M)-1][2]))
file.write('\n')

# Print normal vectors to file
for i in range(len(M)-1, -1, -1):
    file.write(str(M[i][0]))
    file.write(' ')
    file.write(str(M[i][1]))
    file.write(' ')
    file.write(str(M[i][2]))
    file.write('\n')
    
file.close()