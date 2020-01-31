#!/usr/bin/python

from cagd.vec import vec2, vec3
from cagd.polyline import polyline
import copy

class bezier_curve:
    def __init__(self, degree):
        assert (degree >= 0)
        self.degree = degree
        self.control_points = [None for i in range(degree + 1)]
        self.color = "black"

    def set_control_point(self, index, val):
        assert (index >= 0 and index <= self.degree)
        self.control_points[index] = val

    def get_control_point(self, index):
        assert (index >= 0 and index <= self.degree)
        return self.control_points[index]

    #evaluates the curve at t
    def evaluate(self, t):
        return self.__de_casteljeau(t, 1)[0]

    #evaluates tangent at t
    def tangent(self, t):
        last_two_ctrl_pts = self.__de_casteljeau(t, 2)
        a = last_two_ctrl_pts[0]
        b = last_two_ctrl_pts[1]
        return b - a

    #calculates the normal at t
    def normal(self, t):
        pass

    #syntactic sugar so bezier curve can be evaluated as curve(t)
    #instead of curve.evaluate(t)
    def __call__(self, t):
        return self.evaluate(t)

    #calculates the de-casteljeau scheme until the column only has stop elements
    def __de_casteljeau(self, t, stop):
        assert (stop >= 1)
        column = self.control_points
        while len(column) > stop:
            new_column = [None for i in range(len(column) -1)]
            for i in range(len(new_column)):
                new_column[i] = (1 - t) * column[i] + t * column[i + 1]
            column = new_column
        return column

    def get_color(self):
        return self.color

    def set_color(self, color):
        self.color = color

    #calculates the bezier representation of the derivative
    def get_derivative(self):
        pass

    def get_axis_aligned_bounding_box(self):
        min_vec = copy.copy(self.control_points[0])
        max_vec = copy.copy(self.control_points[0])
        for p in self.control_points:
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
        p0 = self(0)
        for i in range(1, num_samples + 1):
            t = i / num_samples
            p1 = self(t)
            scene.draw_line(p0, p1, self.color)
            p0 = p1

    def get_polyline_from_control_points(self):
        pl = polyline()
        for p in self.control_points:
            pl.append_point(p)
        return pl


class bezier_surface:
    #creates a bezier surface of degrees n,m
    #the degree parameter is a tuple (n,m)
    def __init__(self, degree):
        d1, d2 = degree
        assert (d1 >= 0 and d2 >= 0)
        self.degree = degree
        self.control_points = [[None for i in range(d2 + 1)] for j in range(d1 + 1)]
        white = (1,1,1)
        self.color = (white, white, white, white)
        self.curveature = (None, None, None, None)

    def set_control_point(self, index1, index2, val):
        assert (index1 >= 0 and index1 <= self.degree[0])
        assert (index2 >= 0 and index2 <= self.degree[1])
        self.control_points[index1][index2] = val

    def get_control_point(self, index1, index2):
        assert (index1 >= 0 and index1 <= self.degree[0])
        assert (index2 >= 0 and index2 <= self.degree[1])
        return self.control_points[index1][index2]

    def evaluate(self, t1, t2):
        return self.__de_casteljeau(t1, t2, (1, 1))[0][0]

    #sets the colors at the corners
    #c00 is the color at u=v=0, c01 is the color at u=0 v=1, etc
    #a color is a tuple (r,g,b) with values between 0 an 1
    def set_colors(self, c00, c01, c10, c11):
        self.color = (c00, c01, c10, c11)

    #sets the curveature at the corners
    #c00 is the curveature at u=v=0, c01 is the curveature at u=0 v=1, etc
    def set_curveature(self, c00, c01, c10, c11):
        self.curveature = (c00, c01, c10, c11)

    def __call__(self, t):
        t1, t2 = t
        return self.evaluate(t1, t2)

    def __de_casteljeau(self, t1, t2, stop):
        s1, s2 = stop
        d1, d2 = self.degree
        assert (s1 >= 1 and s2 >= 1)
        d1 += 1 #number of control points in each direction
        d2 += 1

        #apply the casteljeau scheme in one direction,
        #ie, reduce dimension from (d1, d2) to (s1, d2)
        column = self.control_points
        while d1 > s1:
            d1 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t1) * column[i][j] + t1 * column[i + 1][j]
            column = new_column

        #apply the casteljeau scheme in the other direction,
        #ie, reduce dimension from (s1, d2) to (s1, s2)
        while d2 > s2:
            d2 -= 1
            new_column = [[None for i in range(d2)] for j in range(d1)]
            for i in range(d1):
                for j in range(d2):
                    new_column[i][j] = (1 - t2) * column[i][j] + t2 * column[i][j + 1]
            column = new_column

        return column

    def normal(self, t1, t2):
        pass

    def get_derivative(self, direction):
        pass

    def subdivide(self, t1, t2):
        b0,b1 = self.__subdivide_u(t1)
        b00,b01 = b0.__subdivide_v(t2)
        b10,b11 = b1.__subdivide_v(t2)
        return [b00, b01, b10, b11]

    def __subdivide_u(self, t):
        du, dv = self.degree
        left = bezier_surface((du, dv))
        right = bezier_surface((du, dv))
        for k in range(du+1):
            pts = self.__de_casteljeau(t, 0, (du-k+1, dv+1))
            left.control_points[k] = pts[0]
            right.control_points[k] = pts[-1]
        return (left, right)

    def __subdivide_v(self, t):
        du, dv = self.degree
        left = bezier_surface((du, dv))
        right = bezier_surface((du, dv))
        for k in range(dv+1):
            pts = self.__de_casteljeau(0, t, (du+1, dv-k+1))
            for i in range(dv+1):
                left.control_points[i][k] = pts[i][0]
                right.control_points[i][k] = pts[i][-1]
        return (left, right)
        


class bezier_patches:
    CURVEATURE_GAUSSIAN = 0
    CURVEATURE_AVERAGE = 1
    CURVEATURE_PRINCIPAL_MAX = 2 #Maximale Hauptkruemmung
    CURVEATURE_PRINCIPAL_MIN = 3 #Minimale Hauptkruemmung
    COLOR_MAP_LINEAR = 4
    COLOR_MAP_CUT = 5
    COLOR_MAP_CLASSIFICATION = 6

    def __init__(self):
        self.patches = []

    def __len__(self):
        return len(self.patches)

    def __getitem__(self, p):
        return self.patches[p]

    def __setitem__(self, i, p):
        self.patches[i] = p

    def __delitem__(self, p):
        del self.patches[p]

    def __iter__(self):
        return iter(self.patches)

    def append(self, p):
        self.patches.append(p)

    #refines patches by subdividing each patch into four new patches
    #there are 4^num times more patches after calling this function
    def refine(self, num):
        for i in range(num):
            new_patches = bezier_patches()
            for p in self:
                new = p.subdivide(0.5, 0.5)
                for n in new:
                    new_patches.append(n)
            self.patches = new_patches

    def visualize_curveature(self, curveature_mode, color_map):
        import numpy as np
        #calculate curveatures at each corner point

        (m, n) = (4, 4)
        for patch in self:
            b = patch.control_points
            # B(bu)
            bu = [[m * (b[i + 1][j] - b[i][j]) for j in range(n)] for i in range(m - 1)]
            # B(bv)
            bv = [[n * (b[i][j + 1] - b[i][j]) for j in range(n - 1)] for i in range(m)]
            # B(buu)
            buu = [[m * (m - 1) * (b[i + 2][j] - 2 * b[i + 1][j] + b[i][j]) for j in range(n)] for i in range(m - 2)]
            # B(buv)
            buv = [[n * m * (b[i + 1][j + 1] - b[i + 1][j] - b[i][j + 1] + b[i][j]) for j in range(n - 1)] for i in
                   range(m - 1)]
            # B(bvv)
            bvv = [[n * (n - 1) * (b[i][j + 2] - 2 * b[i][j + 1] + b[i][j]) for j in range(n - 2)] for i in range(m)]
        for patch in self:
            w = self.calculate_derivative(np.array(patch.control_points), 1)

        E = self.cal_frobenius_inner_product(bu, bu)
        F = self.cal_frobenius_inner_product(bu, bv)
        G = self.cal_frobenius_inner_product(bv, bv)

        N =

        e = self.cal_frobenius_inner_product(N, buu)
        f = self.cal_frobenius_inner_product(N, buv)
        g = self.cal_frobenius_inner_product(N, bvv)


        #set colors according to color map
        # Die Krümmungswerte müssen entsprechend dem Parameter color_map auf Farbwerte abgebildet werden.

        x = [0.] * 100
        result = [vec3(0,0,0)] * 100
        # Die Funktionen f1 : IR → [0, 1], die jedem Krümmungswert einen Wert aus dem Bereich [0, 1]
        # zuordnen; auf das Ergebnis wird die Hilfsfunktion h angewendet
        for i in range(len(x)):
            if (x[i] < 0):
                result[i] = self.cal_color(0)
            elif (1 < x[i]):
                result[i] = self.cal_color(1)
            else:
                result[i] = self.cal_color(x[i])

        # Mit der Funktion f2 wird κmin auf 0, κmax auf 1 abgebildet und die Werte dazwischen
        # werden linear interpoliert:
        x_min = min(x)
        x_max = max(x)
        for i in range(len(x)):
            result[i] = self.cal_color((x[i] - x_min)/(x_max - x_min))

        # Die Funktion f3, die die hyperbolischen Punkte blau,
        # die parabolischen und die Flachpunkte grün und die elliptischen Punkte rot darstellt
        for i in range(len(x)):
            if (x[i] < 0):
                result[i] = self.cal_color(0)
            elif (1 < x[i]):
                result[i] = self.cal_color(0.5)
            else:
                result[i] = self.cal_color(1)

    # Die Hilfsfunktion h: [0, 1] → [0, 1]3, die 0 auf (0, 0, 1) (blau),
    # 1/2 auf (0, 1, 0) (grün) und 1 auf (1, 0, 0) (rot) abbildet und
    # die die Werte dazwischen linear interpoliert:
    def cal_color(self, x):
        assert (0. <= x <= 1.)
        if (0. <= x <= 0.25):
            h_x = vec3(0, 4*x, 1)
        elif (0.25 < x <= 0.5):
            h_x = vec3(0, 1, 2-4*x)
        elif (0.5 < x <= 0.75):
            h_x = vec3(4*x-2, 1, 0)
        elif (0.75 < x <= 1.):
            h_x = vec3(1, 4-4*x, 0)
        return h_x

    def calculate_derivative(self, control_points, n_derivatives):
        import numpy as np
        w = {0: control_points}
        for i in range(n_derivatives):
            n = len(w[i])
            w[i + 1] = np.array([(n - 1) * (w[i][j + 1] - w[i][j]) for j in range(n - 1)])
        return w

    def cal_frobenius_inner_product(self, a, b):
        assert ((len(a) == len(b)) and (len(a[0]) == len(b[0])))
        c = 0
        for i in range(len(a)):
            for j in range(len(a[0])):
                c += a[i][j]*b[i][j]
        return c




    def export_off(self):
        def export_point(p):
            return str(p.x) + " " + str(p.y) + " " + str(p.z)
        def export_colors(c):
            s = ""
            for x in c:
                s += str(x)
                s += " "
            s += "1" #opacity
            return s

        s = "CBEZ333\n"
        for patch in self:
            #coordinates
            for row in patch.control_points:
                for p in row:
                    s += export_point(p)
                    s += "\n"

            #colors
            s += export_colors(patch.color[0])
            s += "\n"
            s += export_colors(patch.color[2])
            s += "\n"
            s += export_colors(patch.color[1])
            s += "\n"
            s += export_colors(patch.color[3])
            s += "\n"
            s += "\n"

        return s
