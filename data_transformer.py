import re

width_in_mm = 14.7
height_in_mm = 9.8
depth_in_mm = 8.575

voxel_size = 0.02449999935925007

org = open('data/cochlea_spiral_points.txt', 'r')
trn = open('data/cochlea_spiral_points_transformed.txt', 'w')

while 1:
    line_in = org.readline()
    if line_in == '':
        break
    else:
        numbers = re.findall(r'-?\d+.\d+', line_in)
        trn.write(str(round(((float(numbers[0])) + width_in_mm/2)/voxel_size)))
        trn.write(' ')
        trn.write(str(round(((float(numbers[1])) + height_in_mm/2)/voxel_size)))
        trn.write(' ')
        trn.write(str(round(((float(numbers[2])) + depth_in_mm/2)/voxel_size)))
        trn.write('\n')

trn.close()
org.close()







