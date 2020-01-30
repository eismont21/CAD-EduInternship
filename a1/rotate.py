#!/usr/bin/python

from cagd.polyline import polyline
from cagd.spline import spline, spline_surface, knots
from cagd.bezier import bezier_surface, bezier_patches
from cagd.vec import vec2, vec3
import cagd.scene_2d as scene_2d
import math
from math import *


#generates a rotational surface by rotating around the z axis
#the input spline is assumed to be on the xz-plane
#num_samples refers to the number of interpolation points in the rotational direction
#returns a spline surface in three dimensions
def generate_rotation_surface(spl, num_samples):
    ss = spline_surface(3)

    # Jeder Kontrollpunkt des Eingabesplines müssen um die z-Achse rotiert werden
    d = spl.control_points
    c = [[vec3(0, 0, 0)] * num_samples for i in range(len(d))]
    for i in range(len(d)):
        for j in range(num_samples):
            k_x = d[i].x * math.cos(2 * math.pi * j / num_samples)
            k_y = d[i].x * math.sin(2 * math.pi * j / num_samples)
            k_z = d[i].y
            c[i][j] = vec3(k_x, k_y, k_z)

    #Für jedes feste i interpolieren wir die Punkte ci0, . . . , ci,R−1 wie in Versuch 2
    b = [[vec3(0, 0, 0)] * (num_samples+3) for i in range(len(d))]
    for i in range(len(d)):
        pts = [vec2(0, 0)] * num_samples
        for j in range(num_samples):
            pts[j] = vec2(c[i][j].x, c[i][j].y)
        s = spline.interpolate_cubic_periodic(pts)
        #Kontrollpunkte bij ∈ R3. (Für festes i liegen die bij in der Ebene z = zi.)
        for k in range(num_samples+3):
            l_x = s.control_points[k].x
            l_y = s.control_points[k].y
            l_z = d[i].y
            b[i][k] = vec3(l_x, l_y, l_z)

    u = spl.knots
    v = s.knots

    ss.knots = (u, v)
    ss.control_points = b

    return ss

pts = [ vec2(0.05,6),
        vec2(0.2,6),
        vec2(1,5),
        vec2(.25,4),
        vec2(1,3.75),
        vec2(.55,3.5),
        vec2(.5,3.4),
        vec2(.5,.6),
        vec2(.55,.5),
        vec2(0.8, 0.2),
        vec2(.95, 0.1),
        vec2(1,0)]

spl = spline.interpolate_cubic(spline.INTERPOLATION_CHORDAL, pts)
#you can activate these lines to view the input spline
#spl.set_color("#0000ff")
#sc = scene_2d.scene()
#sc.set_resolution(900)
#sc.add_element(spl)
#sc.write_image()
#sc.show()

surface = generate_rotation_surface(spl, 12)

bezier_patches = surface.to_bezier_patches()

path = "surfaces.off"
f = open(path, 'w')
f.write(bezier_patches.export_off())
f.close()
