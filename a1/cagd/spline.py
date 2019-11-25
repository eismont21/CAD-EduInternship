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
        #print(a, b , t)
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
        #print("knot_index = ", knot_index)
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

        return d[(len(d)-stop):]
        #pass

    #adjusts the control points such that it represents the same function,
    #but with an added knot
    def insert_knot(self, t):
        print("knot to insert:", t)
        print("knots before insert: ",self.knots.knots[3:-3])
        #self.knots.insert(t)
        #print("knots after insert: ",self.knots.knots[3:-3], end="\n\n")
        index = self.knots.knot_index(t)

        #ctrl_pts =[]
        #a = (t-self.knots[index-1])/(self.knots[index-1+self.degree] - self.knots[index-1])
        #ctrl_pts.append((vec2(1,1) - a*self.control_points[index - 2])+a*self.control_points[index-1])
        #a = (t-self.knots[index-2])/(self.knots[index-2+self.degree] - self.knots[index-2])
        #ctrl_pts.append((vec2(1,1) - a*self.control_points[index - 3])+a*self.control_points[index-2])
        #a = (t-self.knots[index-3])/(self.knots[index-3+self.degree] - self.knots[index-3])
        #ctrl_pts.append((vec2(1,1) - a*self.control_points[index - 4])+a*self.control_points[index-3])
        ctrl_pts = self.de_boor(t, 3)
        control_pts_x = [p.x for p in ctrl_pts]
        control_pts_y = [p.y for p in ctrl_pts]
        print("controll points to add", control_pts_x, control_pts_y, sep="\n", end="\n\n")
        #ctrl_pts = self.de_boor(t, 3)
        #control_pts_x = [p.x for p in ctrl_pts]
        #control_pts_y = [p.y for p in ctrl_pts]
        #print("controll points to add", control_pts_x, control_pts_y, sep="\n", end="\n\n")

        print("index, len ctrl pts to insert = ", index, len(ctrl_pts), end="\n\n")

        control_pts_x = [p.x for p in self.control_points]
        control_pts_y = [p.y for p in self.control_points]
        print("controll points before insert", control_pts_x, control_pts_y, sep="\n")
        self.control_points = self.control_points[:(index-self.degree+1)] + ctrl_pts + self.control_points[(index):]

        self.knots.insert(t)

        print("knots after insert: ",self.knots.knots[3:-3], end="\n\n")
        control_pts_x = [p.x for p in self.control_points]
        control_pts_y = [p.y for p in self.control_points]
        print("controll points after insert", control_pts_x, control_pts_y, sep="\n", end="\n\n")

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
        s = spline(3)
        s.control_points = points
        n = len(points)

        # Parametrisierung
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

            #chordale
            for i in range(0, n-1):
                d[i] = math.sqrt((points[i+1] - points[i]).x ** 2 + (points[i+1] - points[i]).y ** 2)

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
            #print(t)

        #Knot vector
        m = len(t)
        u = [0.] * (m+6) #corresponds t in Assignment, 6 is adding 3 points left and 3 points right
        u[0], u[1], u[2] = float(t[0]), float(t[0]), float(t[0]) #t1 =t2 = t3
        for i in range(m):
            u[i+3] = float(t[i])
        u[m+3], u[m+4], u[m+5] = float(t[m-1]), float(t[m-1]), float(t[m-1])

        knots_t = knots(len(u))
        knots_t.knots = u
        s.knots = knots_t
        #print("u = ", u)

        #creating of the equation Ax=p
        #A is tridiagonal and consists of main_diag, upper_diag and under_diag
        #p consists of elements of type vec2 namely points

        #calculating p
        p = [0] * (n+2)
        p[0] = points[0]
        p[1] = vec2(0.0, 0.0)
        p[n+1] = points[n-1]
        p[n] = vec2(0.0, 0.0)
        p[2:n] = points[1:(n-1)]

        #calculating A
        main_diag = [0.] * (n + 2)
        under_diag = [0.] * (n + 2)
        upper_diag = [0.] * (n + 2)

        for i in range(2, n):
            #print("i = ", i)
            ai = (u[i+2]-u[i])/(u[i+3] - u[i])
            bi = (u[i+2]-u[i+1])/(u[i+3] - u[i+1])
            ci = (u[i+2]-u[i+1])/(u[i+4] - u[i+1])
            #print("ai, bi, ci = ", ai, bi, ci, sep=" ", end="\n")
            under_diag[i] = (1-bi)*(1-ai)
            main_diag[i] = (1-bi)*ai + bi*(1-ci)
            upper_diag[i] = bi*ci

        main_diag[0] = 1.0
        main_diag[n + 1] = 1.0
        main_diag[1] = 1 + (u[4] - u[2]) / (u[5] - u[2])
        main_diag[n] = - (u[n+1] - u[n]) / (u[n+3] - u[n]) + 2

        upper_diag[n] = -1.0
        upper_diag[n+1] = .0
        upper_diag[0] = .0
        upper_diag[1] = - (u[4] - u[2]) / (u[5] - u[2])

        under_diag[0] = .0
        under_diag[1] = -1.0
        under_diag[n+1] = .0
        under_diag[n] = -1.0 + (u[n+1] - u[n]) / (u[n+3] - u[n])

        p_x = [el.x for el in p]
        p_y = [el.y for el in p]
        #print("main upper under p_x p_y:", main_diag, upper_diag, under_diag, p_x, p_y, sep="\n")

        #solve equation Ax = p
        x = utils.solve_tridiagonal_equation(under_diag, main_diag, upper_diag, p) # divion by 0!!!
        s.control_points = x
        return s

    #generates a spline that interpolates the given points and fulfills the definition
    #of a periodic spline
    #returns that spline object
    def interpolate_cubic_periodic(points):
        s = spline(3)
        n = len(points) + 3 # len(points) + degree

        # Knot vector
        u = [0] * (n+4) # m(=n-1) + degree(=3) + 1 + 1
        for i in range(n+4):
            u[i] = i
        knots_t = knots(len(u))
        knots_t.knots = u
        s.knots = knots_t

        # calculating A
        main_diag = [2/3] * (n-3)
        under_diag = [1/6] * (n-3)
        upper_diag = [1/6] * (n-3)

        # solve equation Ax = p
        x = utils.solve_almost_tridiagonal_equation(under_diag, main_diag, upper_diag, points)
        x.insert(0, x[len(x)-1]) # d[0] = d[m]
        x.append(x[1]) # d[m+1] = d[1]
        x.append(x[2]) # d[m+2] = d[2]
        s.control_points = x
        return s

    #for splines of degree 3, generate a parallel spline with distance dist
    #the returned spline is off from the exact parallel by at most eps
    def generate_parallel(self, dist, eps):
        assert(self.degree == 3)
        #s_parallel = spline(3)
        #s_parallel.knots = self.knots
        pts = []

        knotss = [self(t) for t in self.knots]
        for i in range(len(knotss)):
            #print("p = ", p.x, p.y)
            tang = self.tangent(self.knots[i])
            x = (tang.y/sqrt(tang.x**2 + tang.y**2))*dist
            y = -(tang.x/sqrt(tang.x**2 + tang.y**2))*dist
            pts.append(vec2(knotss[i].x + x, knotss[i].y + y))

        pts_x = [p.x for p in pts]
        pts_y = [p.y for p in pts]
        #print(pts_x, pts_y, sep="\n")
        s = spline.interpolate_cubic(self.INTERPOLATION_CHORDAL, pts[3:-3])

        #"Setzen Sie daher die Knoten des parallelen Splines auf die Knoten des Eingabesplines".
        self.knots = s.knots
        #knotss = [self(t) for t in self.knots]
        knotss_new = [s(t) for t in s.knots]
        for i in range(len(knotss_new)-1):
            #old spline
            middle_x = (knotss[i+1].x + knotss[i].x)/2
            middle_y = (knotss[i + 1].y + knotss[i].y) / 2
            #new spline
            middle_x2 = (knotss_new[i + 1].x + knotss_new[i].x) / 2
            middle_y2 = (knotss_new[i+1].y + knotss_new[i].y)/2
            #Jetzt soll in der Mitte zwischen zwei Knoten die Distanz der Splines berechnet werden.
            middle_dist = sqrt((middle_x-middle_x2)**2+(middle_y-middle_y2)**2)

            if (abs(middle_dist-dist) > eps):
                #Bei einer Abweichung von mehr als eps von der geforderten Distanz soll an dieser Stelle
                #ein neuer Knoten in den originalen Spline eingefügt werden um eine bessere Approximation zu gewinnen.

                # README
                # если ты раскоммитишь след строчку добавления, то получишь ошибку питона, я хз почему, уже заебавси
                self.insert_knot((self.knots[i+1] + self.knots[i])/2)
                #continue

                #print(knotss[i].x, knotss_new[i].x, knotss[i].y, knotss_new[i].y, abs(middle_dist-dist))
        knotss = [self(t) for t in self.knots[3:-3]]
        pts_x = [p.x for p in knotss]
        pts_y = [p.y for p in knotss]
        print(pts_x, pts_y, sep="\n")
        print("self knots", self.knots.knots[3:-3])
        pts = []
        for i in range(len(knotss)):
            #print("p = ", p.x, p.y)
            tang = self.tangent(self.knots[3:-3][i])
            #print("tang = ", tang.x, tang.y)
            x = (tang.y/sqrt(tang.x**2 + tang.y**2))*dist
            y = -(tang.x/sqrt(tang.x**2 + tang.y**2))*dist
            print("x, y = ", x, y)
            pts.append(vec2(knotss[i].x + x, knotss[i].y + y))

        pts_x = [p.x for p in pts]
        pts_y = [p.y for p in pts]
        print(pts_x, pts_y, sep="\n")
        s = spline.interpolate_cubic(self.INTERPOLATION_CHORDAL, pts)
        return s

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
