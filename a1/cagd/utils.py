#!/usr/bin/python
from cagd.vec import vec2

#solves the system of linear equations Ax = res
#where A is a tridiagonal matrix with diag2 representing the main diagonal
#diag1 and diag3 represent the lower and upper diagonal respectively
#all four parameters are vectors of size n
#the first element of diag1 and the last element of diag3 are ignored
#therefore diag1[i], diag2[i] and diag3[i] are located on the same row of A
#res is an array of vec2
#in fact the system is solved for x and y coordinate
#returns the array with elements of type vec2
def solve_tridiagonal_equation(diag1, diag2, diag3, res):
    assert(len(diag1) == len(diag2) == len(diag3) == len(res))

    res_x = [el.x for el in res]
    res_y = [el.y for el in res]
    sol_x = solve_tridiagonal_equation_one_coordinate(diag1, diag2, diag3, res_x)
    sol_y = solve_tridiagonal_equation_one_coordinate(diag1, diag2, diag3, res_y)
    return [vec2(sol_x[i], sol_y[i]) for i in range (len(res))]


def solve_tridiagonal_equation_one_coordinate(diag1, diag2, diag3, res):
    n = len(diag2)
    x = [0] * n
    y = [0] * n
    v = [0] * n
    a = diag1
    z = [0] * n
    c = diag3
    b = diag2
    d = res

    for i in range(n):
        z[i] = 1 / (b[i] - a[i] * v[i])
        if i != n - 1:
            v[i + 1] = z[i] * c[i]
        if i == 0:
            y[i] = z[i] * d[i]
        else:
            y[i] = z[i] * (d[i] - a[i] * y[i - 1])

    x[n - 1] = y[n - 1]
    for i in range(n-2, -1, -1):
        x[i] = y[i] - v[i + 1] * x[i + 1]

    return x
#solves the system of linear equations Ax = res
#where A is an almost tridiagonal matrix with diag2 representing the main diagonal
#diag1 and diag3 represent the lower and upper diagonal respectively
#all four parameters are vectors of size n
#the first element of diag1 and the last element of diag3 represent the top right and bottom left elements of A
#diag1[i], diag2[i] and diag3[i] are located on the same row of A
def solve_almost_tridiagonal_equation(diag1, diag2, diag3, res):
    assert(len(diag1) == len(diag2) == len(diag3) == len(res))
    pass

