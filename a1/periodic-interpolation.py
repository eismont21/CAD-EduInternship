#!/usr/bin/python

from cagd.polyline import polyline
from cagd.spline import spline, knots
from cagd.vec import vec2
import cagd.scene_2d as scene_2d
from math import sin,cos,pi

#returns a list of num_samples points that are uniformly distributed on the unit circle
def unit_circle_points(num_samples):
    a = 2*pi/num_samples
    return [vec2(cos(a*i), sin(a*i)) for i in range(num_samples)]

#calculates the deviation between the given spline and a unit circle
def calculate_circle_deviation(spline):
    return None


#interpolate 6 points with a periodic spline to create the number "8"
#pts = [vec2( 0, 2.5), vec2(-1, 1), vec2( 1,-1), vec2( 0,-2.5), vec2(-1,-1), vec2(1,1)]
#s = spline.interpolate_cubic_periodic(pts)
#p = s.get_polyline_from_control_points()
#p.set_color("blue")
sc = scene_2d.scene()
sc.set_resolution(900)
#sc.add_element(s)
#sc.add_element(p)

#generate a spline that approximates the unit circle
n = 8
circle_pts = unit_circle_points(n)
for el in circle_pts:
    print(el.x, el.y)
circle = spline.interpolate_cubic_periodic(circle_pts)
sc.add_element(circle)
#error = calculate_circle_deviation(circle)
#print("The error is: " + str(error))

sc.write_image()
sc.show()
