#!/usr/bin/python

from cagd.spline import spline_surface, knots
from cagd.bezier import bezier_surface, bezier_patches
from cagd.vec import vec3


f = lambda x,y: x*x/5 + x*y/4 + 10/(1+x*x+y*y) + y/2
ctrl_pts = [[vec3(x, y, f(x,y)) for x in range(-5, 5)] for y in range(-5, 5)]

d = 3
m = d + len(ctrl_pts) + 1
ku = knots(m)
kv = knots(m)
for i in range(m):
    ku[i] = i
    kv[i] = i

surface = spline_surface((d,d))
surface.control_points = ctrl_pts
surface.knots = (ku, kv)

bezier_patches = surface.to_bezier_patches()
bezier_patches.refine(1) #refine surface into more patches for more detailed coloring
bezier_patches.visualize_curveature(bezier_patches.CURVEATURE_GAUSSIAN, bezier_patches.COLOR_MAP_LINEAR)

path = "surfaces.off"
f = open(path, 'w')
f.write(bezier_patches.export_off())
f.close()
