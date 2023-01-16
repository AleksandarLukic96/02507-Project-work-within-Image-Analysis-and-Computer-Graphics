import numpy as np

def trilinearInterpolation(x, y, z, img_data):
    dx = x - int(x)
    dy = y - int(y)
    dz = z - int(z)

    v1 = img_data(np.ceil(x), np.ceil(y), np.ceil(z))
    v2 = img_data(np.ceil(x), np.ceil(y), np.floor(z))
    v3 = img_data(np.ceil(x), np.floor(y), np.floor(z))
    v4 = img_data(np.ceil(x), np.floor(y), np.ceil(z))
    v5 = img_data(np.floor(x), np.ceil(y), np.ceil(z))
    v6 = img_data(np.floor(x), np.floor(y), np.ceil(z))
    v7 = img_data(np.floor(x), np.ceil(y), np.floor(z))
    v8 = img_data(np.floor(x), np.floor(y), np.floor(z))

    # Upper square: bilinear interpolation
    g1 = bilinearInterpolation(v1, v2, v3, v4, dx, dy)

    # Lower square: bilinear interpolation
    g2 = bilinearInterpolation(v5, v6, v7, v8, dx, dy)

    # Linear interpolation
    voxelValue = g1 * (1 - dz) \
                 + g2 * (dz)

    return voxelValue


def bilinearInterpolation(v1, v2, v3, v4, dx, dy):
    g = v1 * (1 - dx) * (1 - dy) \
         + v2 * (dx) * (1 - dy) \
         + v3 * (1 - dx) * (dy) \
         + v4 * (dx * dy)
    return g