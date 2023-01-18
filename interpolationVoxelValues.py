import numpy as np

def trilinearInterpolation(point, img_data):
    x = point[0]
    y = point[1]
    z = point[2]

    dx = x - int(x)
    dy = y - int(y)
    dz = z - int(z)

    cx = int(np.ceil(x))
    cy = int(np.ceil(y))
    cz = int(np.ceil(z))
    fx = int(np.floor(x))
    fy = int(np.floor(y))
    fz = int(np.floor(z))
    
    v1 = img_data[cx, cy, cz]
    v2 = img_data[cx, cy, fz]
    v3 = img_data[cx, fy, fz]
    v4 = img_data[cx, fy, cz]
    v5 = img_data[fx, cy, cz]
    v6 = img_data[fx, fy, cz]
    v7 = img_data[fx, cy, fz]
    v8 = img_data[fx, fy, fz]

    # Upper square: bilinear interpolation
    g1 = bilinearInterpolation(v1, v2, v3, v4, dx, dy)

    # Lower square: bilinear interpolation
    g2 = bilinearInterpolation(v5, v6, v7, v8, dx, dy)

    # Linear interpolation
    voxel_value = g1 * (1 - dz) + g2 * dz

    return voxel_value


def bilinearInterpolation(v1, v2, v3, v4, dx, dy):
    g = v1 * (1 - dx) * (1 - dy) \
         + v2 * dx * (1 - dy) \
         + v3 * (1 - dx) * dy \
         + v4 * dx * dy
    return g