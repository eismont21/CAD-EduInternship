#! /usr/bin/python

from cagd.vec import vec2, vec3
from cagd.polyline import polyline
import cagd.utils as utils
import copy
import math
from math import *

class spline:
    #Interpolation modes
    INTERPOLATION_EQUIDISTANT = 0
    INTERPOLATION_CHORDAL = 1
    INTERPOLATION_CENTRIPETAL = 2
    INTERPOLATION_FOLEY = 3

    def __init__(self, degree):
        assert(degree >= 1)
        self.degree = degree
        self.knots = None
        self.control_points = []
        self.color = "black"

    #checks if the number of knots, controlpoints and degree define a valid spline
    def validate(self):
        knots = self.knots.validate()
        points = len(self.knots) == len(self.control_points) + self.degree + 1
        return knots and points

    def evaluate(self, t):
        a, b = self.support()
        assert(a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            #the spline is only defined on the interval [a, b)
            #it is useful to define self(b) as lim t->b self(t)
            t = t - 0.000001
        return self.de_boor(t, 1)[0]

    #returns the interval [a, b) on which the spline is supported
    def support(self):
        return (self.knots[self.degree], self.knots[len(self.knots) - self.degree - 1])

    def __call__(self, t):
        return self.evaluate(t)

    def tangent(self, t):
        a, b = self.support()
        assert(a <= t <= b)
        if t == self.knots[len(self.knots) - self.degree - 1]:
            #the spline is only defined on the interval [a, b)
            #it is useful to define self(b) as lim t->b self(t)
            t = t - 0.000001
        last_two_points = self.de_boor(t, 2)
        return last_two_points[1] - last_two_points[0]

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    #calculates the de_boor scheme at a given value t
    #stops when the column is only "stop" elements long
    #returns that column as a list
    def de_boor(self, t, stop):
        knot_index = self.knots.knot_index(t)
        print("knot_index = ", knot_index)
        d = [self.control_points[i + knot_index - self.degree] for i in range(self.degree + 1)]
        #print(len(d));
        for k in range(1, self.degree+1):
            for j in range(self.degree, k-1, -1):
                alphakj = (t - self.knots[j + knot_index - self.degree]) / (self.knots[j + 1 + knot_index - k] - self.knots[j + knot_index - self.degree])
                d[j] = (1.0 - alphakj) * d[j-1] + alphakj * d[j]
            if self.degree - (k-1) == stop:
                break
        #for el in d:
        #    print(el)
        #    print(" ")

        return d[(len(d)-stop-1):]
        #pass

    #adjusts the control points such that it represents the same function,
    #but with an added knot
    def insert_knot(self, t):
        pass

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
            #print("comparing {0} to {1} and {2}".format(p, min_vec, max_vec))
            if p.x < min_vec.x:
                min_vec.x = p.x
            if p.y < min_vec.y:
                min_vec.y = p.y
            if p.x > max_vec.x:
                max_vec.x = p.x
            if p.y > max_vec.y:
                max_vec.y = p.y
        return (min_vec, max_vec)

    def draw(self, scene, num_samples):
        i = self.degree - 1
        while i < len(self.knots) - self.degree - 2:
            i += 1
            k0 = self.knots[i]
            k1 = self.knots[i+1]
            if k0 == k1:
                continue
            p0 = self(k0)
            for j in range(1, num_samples + 1):
                t = k0 + j / num_samples * (k1 - k0)
                p1 = self(t)
                scene.draw_line(p0, p1, self.color)
                p0 = p1

    def get_polyline_from_control_points(self):
        pl = polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl
            
    #generates a spline that interpolates the given points using the given mode
    #returns that spline object
    def interpolate_cubic(mode, points):
        #Parametrisierung
        s = spline(3)
        s.control_points = points
        n = len(points)
        t = [0.] * n

        if mode == 0:
            for i in range(n):
                t[i] = i

        elif mode == 1:
            for i in range(1, n):
                t[i] = math.sqrt((points[i]-points[i-1]).x ** 2 + (points[i]-points[i-1]).y ** 2) + t[i-1]

        elif mode == 2:
            for i in range(1, n):
                t[i] = ((points[i] - points[i - 1]).x ** 2 + (points[i] - points[i - 1]).y ** 2) ** 0.25 + t[i - 1]

        elif mode == 3:
            d = [0.] * n

            for i in range(1, n):
                d[i] = math.sqrt((points[i] - points[i - 1]).x ** 2 + (points[i] - points[i - 1]).y ** 2) + d[i - 1]

            alpha = [0.] * n

            for i in range(1, n-1):
                a1 = points[i-1].y - points[i].y
                b1 = points[i].x - points[i-1].x
                a2 = points[i].y - points[i+1].y
                b2 = points[i+1].x - points[i].x
                angle = math.acos((b1*b2 + a1*a2)/(math.sqrt(b1**2 + a1**2) * math.sqrt(b2**2 + a2**2)))
                alpha[i] = min(math.pi - angle, math.pi/2)

            for i in range(1, n):
                if i == 1:
                    k = 0
                else:
                    k = 3/2 * alpha[i-1]*d[i-2]/(d[i-2]+d[i-1])
                l = 3/2 * alpha[i]*d[i]/(d[i]+d[i-1])
                t[i] = d[i-1] * (1 + k + l) + t[i-1]

        #Knot vector
        u = [0.] * (n+5)
        for i in range(n):
            u[i+4] = float(t[i]) #индексы хз, в (8) они с 0 начинают или как, просто тогда первые четыре всегда 0 будут, потому что т0 всегда 0, сложнаааа
        u[n+4] = float(t[n-1])
        u[n+3] = float(t[n-1])
        u[n+2] = float(t[n-1])
        u[n+1] = float(t[n-1])


        knots_t = knots(len(u)-1)
        knots_t.knots = u[1:]
        s.knots = knots_t
        print("knots = ", u)

        #creating of the matrix
        p = [0.] * (n+1)
        p[0] = s.de_boor(u[3], 1)[0].y
        p[1] = 0
        p[n] = s.de_boor(u[n+3], 1)[0].y
        p[n-1] = 0
        for i in range(2, n):
            print(i, sep=" ")
            p[i-1] = s.de_boor(u[i+2], 1)[0].y
        main_diag = [0.] * (n + 1)
        under_diag = [0.] * (n + 1)
        upper_diag = [0.] * (n + 1)

        main_diag[0] = 1
        main_diag[n] = 1
        upper_diag[n-1] = -1
        under_diag[0] = -1
        for i in range(2, n):
            print("i = ", i)
            ai = (u[i+2]-u[i])/(u[i+3] - u[i])
            bi = (u[i+2]-u[i+1])/(u[i+3] - u[i+1])
            ci = (u[i+2]-u[i+1])/(u[i+4] - u[i+1])
            print("ai", ai, bi, ci, sep="\n")
            under_diag[i] = (1-bi)*(1-ai)
            main_diag[i] = (1-bi)*ai + bi*(1-ci)
            upper_diag[i] = bi*ci

        print("Diags main upper under p:", main_diag, upper_diag, under_diag, p, sep="\n")
        x = utils.solve_tridiagonal_equation(under_diag, main_diag, upper_diag, p)
        return s






    #generates a spline that interpolates the given points and fulfills the definition
    #of a periodic spline
    #returns that spline object
    def interpolate_cubic_periodic(points):
        pass

    #for splines of degree 3, generate a parallel spline with distance dist
    #the returned spline is off from the exact parallel by at most eps
    def generate_parallel(self, dist, eps):
        assert(self.degree == 3)
        pass


class spline_surface:
    #the two directions of the parameter space
    DIR_U = 0
    DIR_V = 1

    #creates a spline of degrees n,m
    #degree is a tuple (n,m)
    def __init__(self, degree):
        du, dv = degree
        assert(du >= 1 and dv >= 1)
        self.degree = degree
        self.knots = (None, None)  #tuple of both knot vectors
        self.control_points = [[]] #2dim array of control points

    #checks if the number of knots, controlpoints and degree define a valid spline
    def validate(self):
        if len(self.control_points) == 0:
            return False
        k1, k2 = self.knots
        d1, d2 = self.degree
        knots = k1.validate() and k2.validate()
        p1 = len(self.control_points)
        p2 = len(self.control_points[0])
        points1 = len(k1) == p1 + d1 + 1
        points2 = len(k2) == p2 + d2 + 1
        return knots and points1 and points2

    def evaluate(self, u, v):
        s1, s2 = self.support()
        a, b = s1
        c, d = s2
        assert(a <= u <= b and c <= v <= v)
        if u == b:
            u = u - 0.000001
        if v == d:
            v = v - 0.000001
        t = (u, v)
        return self.de_boor(t, (1,1))[0][0]

    #return nested tuple ((a,b), (c,d))
    #the spline is supported in (u,v) \in [a,b)x[c,d]
    def support(self):
        k1, k2 = self.knots
        d1, d2 = self.degree
        s1 = (k1[d1], k1[len(k1) - d1 - 1])
        s2 = (k2[d2], k2[len(k2) - d2 - 1])
        return (s1, s2)

    def __call__(self, u, v):
        return self.evaluate(u, v)

    #calculates the de boor scheme at t = (u,v)
    #until there are only stop = (s1, s2) elements left
    def de_boor(self, t, stop):
        d1, d2 = self.degree
        k1, k2 = self.knots
        s1, s2 = stop
        u, v = t
        m1 = len(self.control_points)
        m2 = len(self.control_points[0])
        
        new_rows = [None for i in range(m1)]
        for row in range(m1):
            spl = spline(d2)
            spl.knots = k2
            spl.control_points = self.control_points[row]
            new_rows[row] = spl.de_boor(v, s2)

        new_pts = [None for i in range(s2)]
        for col in range(s2):
            spl = spline(d1)
            spl.knots = k1
            ctrl_pts = [new_rows[i][col] for i in range(m1)]
            spl.control_points = ctrl_pts
            new_pts[col] = spl.de_boor(u, s1)

        return new_pts

    def insert_knot(self, t):
        pass

class knots:
    #creates a knots array with n elements
    def __init__(self, n):
        self.knots = [None for i in range(n)]

    def validate(self):
        prev = None
        for k in self.knots:
            if k is None:
                return False
            if prev is None:
                prev = k
            else:
                if k < prev:
                    return False
        return True 

    def __len__(self):
        return len(self.knots)

    def __getitem__(self, i):
        return self.knots[i]

    def __setitem__(self, i, v):
        self.knots[i] = v

    def __delitem__(self, i):
        del self.knots[i]

    def __iter__(self):
        return iter(self.knots)

    def insert(self, t):
        i = 0
        while self[i] < t:
            i += 1
        self.knots.insert(i, t)

    def knot_index(self, v):
        n = len(self.knots)
        for i in range(n-1):
            if self.knots[i] <= v < self.knots[i + 1]:
                return i
        return self.knots.index(max(self.knots))-1
        #pass
