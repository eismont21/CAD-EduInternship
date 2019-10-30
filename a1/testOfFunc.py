from cagd.polyline import polyline
from cagd.spline import spline, knots
from cagd.vec import vec2
import cagd.scene_2d as scene_2d

s = spline(3)
s.control_points = [vec2(1,0), vec2(0,3), vec2(2,6), vec2(5,8), vec2(9,6), vec2(10,1), vec2(6,0)]
s.knots = knots(10)
s.knots.knots = [0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0]
for i in range(4, 0, -1):
    print("The column:", i)
    for el in s.de_boor(0.4, i):
        print("x = ", el.x, "; y = ", el.y, end="\n")
