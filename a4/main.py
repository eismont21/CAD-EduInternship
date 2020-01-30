from cube import CubeEdges, CubeVertices, CubeEdgeFlags, CubeTriangles
from a1.cagd.vec import vec3
import math

'''The class representing traingles'''


class Triangle:
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def map(self, f):
        return Triangle(f(self.v1), f(self.v2), f(self.v3))


"""The class representing vertices and faces between those vertices."""


class Mesh:
    def __init__(self, verts=None, faces=None):
        self.verts = verts or []
        self.faces = faces or []

    def extend(self, other):
        l = len(self.verts)
        f = lambda v: v + l
        self.verts.extend(other.verts)
        self.faces.extend(face.map(f) for face in other.faces)

    def export_in_off(self, fileName):
        with open(fileName + ".off", "w") as f:
            f.write("OFF\n")
            f.write("{} {} {}\n".format(len(mesh.verts), len(mesh.faces), 0))
            for v in mesh.verts:
                f.write("{} {} {}\n".format(v.x, v.y, v.z))
            for face in mesh.faces:
                f.write("3 {} {} {}\n".format(face.v1 - 1, face.v2 - 1, face.v3 - 1))

    def export_in_obj(self, fileName):
        with open(fileName + ".obj", "w") as f:
            for v in mesh.verts:
                f.write("v {} {} {}\n".format(v.x, v.y, v.z))
            for face in mesh.faces:
                f.write("f {} {} {}\n".format(face.v1, face.v2, face.v3))


# Default bounds to evaluate
XMIN = -3
XMAX = 3
YMIN = -3
YMAX = 3
ZMIN = -3
ZMAX = 3
ISOVALUE = 2.5


def edge_to_boundary_vertex(edge, f_eval, xyz):
    """Returns the vertex in the middle of the specified edge"""
    p1, p2 = CubeEdges[edge]
    v1 = f_eval[p1]
    v2 = f_eval[p2]
    p1 = vec3(CubeVertices[p1][0], CubeVertices[p1][1], CubeVertices[p1][2])
    p2 = vec3(CubeVertices[p2][0], CubeVertices[p2][1], CubeVertices[p2][2])
    q = (v1*p2 - v2*p1)/(v1 - v2)
    return xyz + q


def marching_cubes_3d_single_cell(f, xyz, k):
    # Evaluate f on each vertex of the cube
    f_eval = [f(xyz + vec3(v[0], v[1], v[2])) for v in CubeVertices]
    # Determine which case it is
    case = sum(2 ** v for v in range(8) if f_eval[v] > 0)
    # what faces are needed (in terms of edges)
    faces = CubeTriangles[case]

    output_verts = []
    output_tris = []

    for i in range(0, len(faces), 3):
        if faces[i] == -1:
            break
        edges = [faces[i:][j] for j in range(3)]
        verts = [edge_to_boundary_vertex(edge, f_eval, xyz) for edge in edges]

        next_vert_index = len(output_verts) + 1
        tri = Triangle(next_vert_index, next_vert_index + 1, next_vert_index + 2)
        output_verts.extend(verts)
        output_tris.append(tri)

    return Mesh(output_verts, output_tris)


def marching_cubes_3d(f):
    """Iterates over a cells of size one between the specified range, and evaluates f to produce
        a boundary by Marching Cubes. Returns a Mesh object."""
    # For each cube, evaluate independently.
    mesh = Mesh()
    #step 0.1
    k = 10
    x_dec_interval = [x / k for x in range(XMIN*k, XMAX*k)]
    y_dec_interval = [x / k for x in range(YMIN*k, YMAX*k)]
    z_dec_interval = [x / k for x in range(ZMIN*k, ZMAX*k)]

    for x in x_dec_interval:
        for y in y_dec_interval:
            for z in z_dec_interval:
                cell_mesh = marching_cubes_3d_single_cell(f, vec3(x, y, z), k)
                mesh.extend(cell_mesh)
    return mesh


def circle(v):
    return math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z) - ISOVALUE


def octaeder(v):
    return abs(v.x) + abs(v.y) + abs(v.z) - ISOVALUE


def torus(v):
    return (v.x ** 2 + v.y ** 2 + v.z ** 2 + 0.375) ** 2 - 2 * (v.x ** 2 + v.y ** 2)


def cube(v):
    return max(abs(v.x), abs(v.y), abs(v.z)) - ISOVALUE


if __name__ == "__main__":
    mesh = marching_cubes_3d(circle)
    fileName = "circle"
    mesh.export_in_off(fileName)
    mesh.export_in_obj(fileName)
