import nibabel as nib
from matplotlib import pyplot as plt
import os

# Define path to directory on any machine
path_to_nii = os.getcwd() + "\\02507\\data\\cochlea\\"
#print("path_to_nii: " + path_to_nii)

# Choose .nii file for examination 
nib_file = "cochlea.nii"

#
epi_img = nib.load(path_to_nii + nib_file)
epi_img_data = epi_img.get_fdata()
epi_img_data.shape

# Function to display slices from the 3 axis
def show_slices(slices):
   """ Function to display row of image slices """
   fig, axes = plt.subplots(1, len(slices))
   for i, slice in enumerate(slices):
       axes[i].imshow(slice.T, cmap = "gray", origin = "lower")

# Set voxel coordinates
slice_x = 100
slice_y = 100
slice_z = 100

# Create slices from voxel coordinates
slice_0 = epi_img_data[slice_x, :, :]
slice_1 = epi_img_data[:, slice_y, :]
slice_2 = epi_img_data[:, :, slice_z]
show_slices([slice_0, slice_1, slice_2])
plt.suptitle("Center slices for EPI image") 

# Save slices to .png file in directory
plt.savefig(os.getcwd() + '\\scliced_img.png')
