# Nice Links
# https://easings.net/

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

def rgb_to_hex(rgb):
    r,g,b = rgb
    return '#%02x%02x%02x' % (r,g,b)

class ArrowScene(Scene):
    def construct(self):
        arrow_left = Arrow(start = RIGHT, end = LEFT, color = GREEN, fill_opacity = 1.0).shift(LEFT)
        arrow_right = Arrow(start = LEFT, end = RIGHT, color = RED, fill_opacity = 1.0).shift(RIGHT)
        arrow_up = Arrow(start = DOWN, end = UP, color = BLUE, fill_opacity = 1.0).shift(UP)
        arrow_down = Arrow(start = UP, end = DOWN, color = YELLOW, fill_opacity = 1.0).shift(DOWN)
        self.next_section("Render arrow LEFT")
        self.play(Create(arrow_left), run_time = 1)
        #self.play(
        #    arrow_left.animate.set_fill(GREEN,opacity=0.5)
        #)
        #self.wait()
        self.next_section("Render arrow UP")
        self.play(Create(arrow_up), run_time = 0.8)
        #self.play(
        #    arrow_up.animate.set_fill(BLUE,opacity=0.5)
        #)
        #self.wait()
        self.next_section("Render arrow RIGHT")
        self.play(Create(arrow_right), run_time = 0.5)
        #self.play(
        #    arrow_right.animate.set_fill(RED,opacity=0.5)
        #)
        #self.wait()
        self.next_section("Render arrow DOWN")
        self.play(Create(arrow_down), run_time = 0.2)
        #self.play(
        #    arrow_down.animate.set_fill(YELLOW,opacity=0.5)
        #)
        self.wait(1)
        self.next_section("Render move arrows")
        self.play(
            arrow_left.animate.shift(UP),
            arrow_up.animate.shift(RIGHT),
            arrow_right.animate.shift(DOWN),
            arrow_down.animate.shift(LEFT)
        )
        self.wait(1)
        self.next_section("Render rotation of arrows")
        self.play(
            arrow_left.animate.shift(0).rotate(-PI / 4),
            arrow_up.animate.shift(0).rotate(-PI / 4),
            arrow_right.animate.shift(0).rotate(-PI / 4),
            arrow_down.animate.shift(0).rotate(-PI / 4)
        )
        self.wait(1)
        self.next_section("Render arrows back to original position and color")
        self.play(
            arrow_left.animate.shift(RIGHT).rotate(-PI / 4).set_color(BLUE),
            arrow_up.animate.shift(DOWN).rotate(-PI / 4).set_color(RED),
            arrow_right.animate.shift(LEFT).rotate(-PI / 4).set_color(YELLOW),
            arrow_down.animate.shift(UP).rotate(-PI / 4).set_color(GREEN)
        )
        self.wait(1)
        self.next_section("Render rotation around center point")
        self.play(
            *[Rotate(mob, angle=2*PI, about_point = ORIGIN, rate_func = rate_functions.ease_in_out_quart, run_time = 2.0) for mob in self.mobjects]
        )
        self.wait(1)
        self.next_section("Render fade out")
        self.play(
            *[FadeOut(mob) for mob in self.mobjects]
            #FadeOut(arrow_left)
            #FadeOut(arrow_up)
            #FadeOut(arrow_right)
            #FadeOut(arrow_down)
        )

class SpiralPointsAndPlanes(ThreeDScene):
    def construct(self):
        # Default camera settings/angles
        phi_init = 65 * DEGREES
        theta_init = -45 * DEGREES

        ################################################################ TODO: Implement a way to render the slices onto Surface objects 
        num = -1
        if num >= 0:
            # Slice to parce as texture onto plane surfaces
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
        colors = [PURE_BLUE, PURE_RED]
        all_colors = color_gradient(colors, data_points_len)
        
        # Make array with spiral points
        points = []
        for i in range(0, data_points_len, 1):
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
        
        index = 55
        x0 = x[index]
        y0 = y[index]
        z0 = z[index]
        i1 = i_norm_vecs[index][0]
        i2 = i_norm_vecs[index][1]
        i3 = i_norm_vecs[index][2]
        j1 = j_norm_vecs[index][0]
        j2 = j_norm_vecs[index][1]
        j3 = j_norm_vecs[index][2]
        
        # Create object
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
        sys.stdout.write("\rCreating slice Completed!\n")
        
        # Rendering:
        sys.stdout.write("\rRendering starting:\n")
        rate = 3
        self.add(center_point)
        self.play(
            *[Create(p) for p in points]
        )
        self.add(plane_surface1)
        self.begin_ambient_camera_rotation(rate = PI/rate)
        self.wait(rate * 2)
        self.stop_ambient_camera_rotation()
        sys.stdout.write("\rRendering Completed!\n")        
        
class SpiralPoints3DExample(ThreeDScene):
    def construct(self):
        phi_init = 65 * DEGREES
        theta_init = -45 * DEGREES

        # Data points to plot
        data = np.genfromtxt(path_spiral_point, dtype = 'float')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]

        # Calculate mean center coordinates for spiral points
        center_x = round(np.mean(x))
        center_y = round(np.mean(y))
        center_z = round(np.mean(z))
        
        # Defines limits of axes according to mean center
        axes = ThreeDAxes(
            x_range = (center_x - 200, center_x + 200, 10),
            y_range = (center_y - 200, center_y + 200, 10),
            z_range = (center_z - 200, center_z + 200, 10)
        )
        
        center_point = Dot3D(point = axes.coords_to_point(center_x, center_y, center_z), color = PURE_GREEN)

        # Set camera orientation and frame center to the center point
        self.set_camera_orientation(phi = phi_init, theta = theta_init, frame_center = center_point)
                
        # Define colors for gradient of rendered points
        colors = [PURE_BLUE, PURE_RED]
        all_colors = color_gradient(colors, len(data))
        
        # Make array with spiral points
        points = []
        for i in range(0, len(data), 1):
        #for i in range(0, len(data), 10):
            point_spiral = Dot3D(point = axes.coords_to_point(x[i], y[i], z[i]), color = all_colors[i])
            points.append(point_spiral)
        
        # Rendering:
        self.begin_ambient_camera_rotation(rate = PI/4)
        
        self.add(center_point)
        self.play(
            *[FadeIn(p) for p in points]
        )

        self.wait(8)
        self.stop_ambient_camera_rotation()

with tempconfig({"quality": "low_quality", "preview": True}):
    scene = SpiralPointsAndPlanes()
    #scene = SpiralPoints3DExample()
    scene.render()