# Nice Links
# https://easings.net/
# https://convertio.co/download/010d0d358fd4073a18e77bef156915533bd43c/
# https://ezgif.com/split/ezgif-2-eb1280f8c1.gif

# Run in termial:
#> manim -p render_3d.py <name of scene>
# Additionally add flags -pqh for high quality

import os
import sys
from manim import *
import numpy as np
import itertools as it
import nibabel as nib
from matplotlib import pyplot as plt
from skimage import img_as_ubyte
from skimage.color import gray2rgb

path = os.getcwd()
path_spiral_point = path + "\\data\\cochlea_spiral_points_transformed.txt"
path_interpolated_points = path + "\\interpolatedPoints.txt"
path_normal_vectors = path + "\\normalVec.txt"
path_saved_slices = path + "\\data\\saved_slices\\"

path_img_test = path + "\\data\\img_test.png"
path_cochlea_tunnel = path + "\\data\\cochlea_tunnel.nii"

# Formatting rgb int values to hexcode as string 
def rgb_to_hex(rgb):
    r,g,b = rgb
    return '#%02x%02x%02x' % (r,g,b)

class SpiralPointsAndPlanes(ThreeDScene):
    def construct(self):
        # Default camera settings/angles
        phi_init = 65 * DEGREES
        theta_init = -45 * DEGREES

        ################################################################ TODO: Implement a functioning way to render the slices onto Surface objects 
        num = -1
        if num >= 0:
            # Slice to parse as texture onto plane surfaces
            slice_0 = path_saved_slices + "slice0.png"
            
            # Convert png data to list of hex-code colors
            img = plt.imread(slice_0)
            img_rgb = img_as_ubyte(img)
            img_hex = []
            step = 1
            for x in range(0, img_rgb.shape[0], step):
                for y in range(0, img_rgb.shape[1], step):
                    img_hex.append(rgb_to_hex(img_rgb[x][y]))
        ################################################################ TODO:
        
        # Span of points extracted from  
        span = 10
        
        # Data points to plot
        data_points = np.genfromtxt(path_interpolated_points, dtype = 'float')
        x = np.ndarray.round(data_points[:, 0][0::span])
        y = np.ndarray.round(data_points[:, 1][0::span])
        z = np.ndarray.round(data_points[:, 2][0::span])
        
        # Normal vectors to define parametric equations
        data_norm_vecs = np.genfromtxt(path_normal_vectors, dtype = 'float')
        n1 = data_norm_vecs[:, 0][0::span]
        n2 = data_norm_vecs[:, 1][0::span]
        n3 = data_norm_vecs[:, 2][0::span]
        
        # Calculate mean center coordinates for spiral points
        center_x = round(np.mean(x))
        center_y = round(np.mean(y))
        center_z = round(np.mean(z))

        # Amount of points to be plotted
        data_points_len = int(round(len(data_points)/span))
        
        # Defines limits of axes according to mean center
        axes_span = 200
        axes = ThreeDAxes(
            x_range = (center_x - axes_span, center_x + axes_span, span),
            y_range = (center_y - axes_span, center_y + axes_span, span),
            z_range = (center_z - axes_span, center_z + axes_span, span)
        )
        
        # 3d dot of Center point for axes-referrence 
        center_point = Dot3D(point = axes.coords_to_point(center_x, center_y, center_z), color = PURE_GREEN)

        # Set camera orientation and frame center to the center point
        self.set_camera_orientation(phi = phi_init, theta = theta_init, frame_center = center_point)

        # Define colors for gradient of rendered points
        #colors = [PURE_BLUE, PURE_RED]
        colors = [BLUE_E, PURE_BLUE]
        all_colors = color_gradient(colors, data_points_len)
        
        # Make array with spiral points
        points = []
        span_points = 1
        for i in range(0, data_points_len, span_points):
            point_spiral = Dot3D(point = axes.coords_to_point(x[i], y[i], z[i]), color = all_colors[i])
            points.append(point_spiral)
            sys.stdout.write("\rGenerating 3D point no. %i" % i)
            sys.stdout.flush()
        sys.stdout.write("\rGenerating 3D point Copmleted!\n")
        
        # Extract equation variables accordingly with normal vectors and center point
        i_norm_vecs = []
        j_norm_vecs = []
        for i in range(0, len(x)):
            norm_vec = np.array([n1[i], n2[i], n3[i]])
            i_vec = np.cross(norm_vec, np.array([0, 0, 1]))
            if all(v == 0 for v in i_vec):
                i_vec = norm_vec
            i_vec_norm = i_vec / np.linalg.norm(i_vec)
            j_vec = np.cross(norm_vec, i_vec)
            j_vec_norm = j_vec / np.linalg.norm(j_vec)
            i_norm_vecs.append(i_vec_norm)
            j_norm_vecs.append(j_vec_norm)
            sys.stdout.write("\rCalculating ortogonal normal vector pair no. %i" % i)
            sys.stdout.flush()
        sys.stdout.write("\rCalculating ortogonal normal vector pair Completed!\n")
        
        create_all_slices = 1
        create_slice1 = 0
        if create_slice1 > 0 or create_all_slices > 0:
            # Calculating parameters for plane1
            index = 22
            x0 = x[index]
            y0 = y[index]
            z0 = z[index]
            i1 = i_norm_vecs[index][0]
            i2 = i_norm_vecs[index][1]
            i3 = i_norm_vecs[index][2]
            j1 = j_norm_vecs[index][0]
            j2 = j_norm_vecs[index][1]
            j3 = j_norm_vecs[index][2]
        
            # Create plane for slice
            sys.stdout.write("\rCreating slice at index %i" % index)
            plane_surface1 = Surface(
                lambda u, v: axes.c2p(
                    x0 + (u * i1) + (v * j1), 
                    y0 + (u * i2) + (v * j2), 
                    z0 + (u * i3) + (v * j3)
                ),
                u_range = [-25, 25],
                v_range = [-25, 25],
                fill_color = WHITE
            )
            sys.stdout.flush()
            sys.stdout.write("\rCreating slice %i Completed!\n" % index)
        
        create_slice2 = 0
        if create_slice2 > 0 or create_all_slices > 0:
            # Calculating parameters for plane1
            index = 44
            x0 = x[index]
            y0 = y[index]
            z0 = z[index]
            i1 = i_norm_vecs[index][0]
            i2 = i_norm_vecs[index][1]
            i3 = i_norm_vecs[index][2]
            j1 = j_norm_vecs[index][0]
            j2 = j_norm_vecs[index][1]
            j3 = j_norm_vecs[index][2]
        
            # Create plane for slice
            sys.stdout.write("\rCreating slice at index %i" % index)
            plane_surface2 = Surface(
                lambda u, v: axes.c2p(
                    x0 + (u * i1) + (v * j1), 
                    y0 + (u * i2) + (v * j2), 
                    z0 + (u * i3) + (v * j3)
                ),
                u_range = [-25, 25],
                v_range = [-25, 25],
                fill_color = WHITE
            )
            sys.stdout.flush()
            sys.stdout.write("\rCreating slice %i Completed!\n" % index)

        create_slice3 = 0
        if create_slice3 > 0 or create_all_slices > 0:
            # Calculating parameters for plane1
            index = 66
            x0 = x[index]
            y0 = y[index]
            z0 = z[index]
            i1 = i_norm_vecs[index][0]
            i2 = i_norm_vecs[index][1]
            i3 = i_norm_vecs[index][2]
            j1 = j_norm_vecs[index][0]
            j2 = j_norm_vecs[index][1]
            j3 = j_norm_vecs[index][2]
        
            # Create plane for slice
            sys.stdout.write("\rCreating slice at index %i" % index)
            plane_surface3 = Surface(
                lambda u, v: axes.c2p(
                    x0 + (u * i1) + (v * j1), 
                    y0 + (u * i2) + (v * j2), 
                    z0 + (u * i3) + (v * j3)
                ),
                u_range = [-25, 25],
                v_range = [-25, 25],
                fill_color = WHITE
            )
            sys.stdout.flush()
            sys.stdout.write("\rCreating slice %i Completed!\n" % index)
        
        create_slice4 = 0
        if create_slice4 > 0 or create_all_slices > 0:
            # Calculating parameters for plane1
            index = 88
            x0 = x[index]
            y0 = y[index]
            z0 = z[index]
            i1 = i_norm_vecs[index][0]
            i2 = i_norm_vecs[index][1]
            i3 = i_norm_vecs[index][2]
            j1 = j_norm_vecs[index][0]
            j2 = j_norm_vecs[index][1]
            j3 = j_norm_vecs[index][2]
        
            # Create plane for slice
            sys.stdout.write("\rCreating slice at index %i" % index)
            plane_surface4 = Surface(
                lambda u, v: axes.c2p(
                    x0 + (u * i1) + (v * j1), 
                    y0 + (u * i2) + (v * j2), 
                    z0 + (u * i3) + (v * j3)
                ),
                u_range = [-25, 25],
                v_range = [-25, 25],
                fill_color = WHITE
            )
            sys.stdout.flush()
            sys.stdout.write("\rCreating slice %i Completed!\n" % index)

        # Rendering:
        sys.stdout.write("\rRendering starting:\n")
        rate = 4
        #self.add(center_point)
        for i in range(0, len(points)):
            sys.stdout.write("\rAdding point[%i]" %i)
            self.add(points[i])    
        sys.stdout.flush()
        sys.stdout.write("\rAdded all points to render!\n")
        self.add(plane_surface1)
        self.add(plane_surface2)
        self.add(plane_surface3)
        self.add(plane_surface4)
        self.begin_ambient_camera_rotation(rate = PI/rate)
        self.wait(rate * 2)
        self.stop_ambient_camera_rotation()
        sys.stdout.write("\rRendering Completed!\n")

with tempconfig({"quality": "low_quality", "preview": True}):
    scene = SpiralPointsAndPlanes()
    #scene = SpiralPoints3DExample()
    scene.render()