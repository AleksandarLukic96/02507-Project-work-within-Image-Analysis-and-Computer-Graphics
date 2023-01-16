# Nice Links
# https://easings.net/

# Run in termial:
#> manim -p render_3d.py <name of scene>
# Additionally add flags -pqh for high quality

import os
from manim import *
import numpy as np

path = os.getcwd()
path_spiral_point = path + "\\data\\cochlea_spiral_points_transformed.txt"

class CreateCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(GREEN, opacity=0.5)  # set the color and transparency
        self.play(Create(circle))  # show the circle on screen
        
class SquareToCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set color and transparency

        square = Square()  # create a square
        square.rotate(PI / 4)  # rotate a certain amount

        self.play(Create(square))  # animate the creation of the square
        self.play(Transform(square, circle))  # interpolate the square into the circle
        self.play(FadeOut(square))  # fade out animation        

class SquareAndCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency

        square = Square()  # create a square
        square.set_fill(BLUE, opacity=0.5)  # set the color and transparency

        square.next_to(circle, RIGHT, buff=0.5)  # set the position
        self.play(Create(circle), Create(square))  # show the shapes on screen

class AnimatedSquareToCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        square = Square()  # create a square

        self.play(Create(square))  # show the square on screen
        self.play(square.animate.rotate(PI / 4))  # rotate the square
        self.play(
            ReplacementTransform(square, circle)
        )  # transform the square into a circle
        self.play(
            circle.animate.set_fill(PINK, opacity=0.5)
        )  # color the circle on screen
        
class DifferentRotations(Scene):
    def construct(self):
        left_square = Square(color=BLUE, fill_opacity=0.7).shift(2 * LEFT)
        right_square = Square(color=GREEN, fill_opacity=0.7).shift(2 * RIGHT)
        self.play(
            left_square.animate.rotate(PI), Rotate(right_square, angle=PI), run_time=2
        )
        self.wait()

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

class ImageImport(Scene):
    def construct(self):
        img = ImageMobject(path + "\\data\\img_test.png")
        self.play(FadeIn(img))
        self.wait(1)
        self.play(FadeOut(img)) 

class Dot3DExample(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(phi = 45*DEGREES, theta = -45*DEGREES, distance = 1000)

        # Data points to plot
        data = np.genfromtxt(path_spiral_point, dtype = 'float')
        x = data[:, 0]
        y = data[:, 1]
        z = data[:, 2]

        axes = ThreeDAxes(
            x_range = (-100, 600, 10),
            y_range = (-100, 600, 10),
            z_range = (-100, 600, 10)
        )
        #for i in range(0, len(data)):
        for i in range(0, 1):
            point_spiral = Dot3D(point=axes.coords_to_point(x[i], y[i], z[i]), color = "#0000FF")
            self.add(axes, point_spiral)

#print(str(len(data)))
#print("x: " + str(x[111]) + " y: " + str(y[0]) + " z: " + str(z[0]))


with tempconfig({"quality": "low_quality", "preview": True}):
    #scene = ArrowScene()
    scene = Dot3DExample()
    scene.render()