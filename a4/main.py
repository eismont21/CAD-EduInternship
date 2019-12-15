"""Provides a function for performing 3D Marching Cubes"""

from cube import CubeEdges, CubeVertices, CubeEdgeFlags, CubeTriangles
from a1.cagd.vec import vec3
import math

class Triangle:
    """A 3d triangle"""
    def __init__(self, v1, v2, v3):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

    def map(self, f):
        return Triangle(f(self.v1), f(self.v2), f(self.v3))

class Mesh:
    """A collection of vertices, and faces between those vertices."""
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
                f.write("3 {} {} {}\n".format(face.v1-1, face.v2-1, face.v3-1))
    def export_in_obj(self, fileName):
        with open(fileName + ".obj", "w") as f:
            for v in mesh.verts:
                f.write("v {} {} {}\n".format(v.x, v.y, v.z))
            for face in mesh.faces:
                f.write("f {} {} {}\n".format(face.v1, face.v2, face.v3))

# Default bounds to evaluate over
XMIN = -3
XMAX = 3
YMIN = -3
YMAX = 3
ZMIN = -3
ZMAX = 3
ISOVALUE= 2.5

def edge_to_boundary_vertex(edge, f_eval, x, y, z):
        """Returns the vertex in the middle of the specified edge"""
        # Find the two vertices specified by this edge, and interpolate between
        # them according to adapt, as in the 2d case
        v0, v1 = CubeEdges[edge]
        f0 = f_eval[v0]
        f1 = f_eval[v1]
        t0 = 1.0 - (0 - f0) / (f1 - f0)
        t1 = 1 - t0
        vert_pos0 = CubeVertices[v0]
        vert_pos1 = CubeVertices[v1]
        return vec3(x + vert_pos0[0] * t0 + vert_pos1[0] * t1,
                  y + vert_pos0[1] * t0 + vert_pos1[1] * t1,
                  z + vert_pos0[2] * t0 + vert_pos1[2] * t1)


def marching_cubes_3d_single_cell(f, x, y, z):
    # Evaluate f on each vertex of the cube
    f_eval = [None] * 8
    for v in range(8):
        v_pos = CubeVertices[v]
        f_eval[v] = f(x + v_pos[0], y + v_pos[1], z + v_pos[2])
    # Determine which case we are
    case = sum(2**v for v in range(8) if f_eval[v] > 0)
    # Ok, what faces do we need (in terms of edges)
    #faces = cases[case]
    faces = CubeTriangles[case]

    output_verts = []
    output_tris = []

    for i in range(0,len(faces),3):
        if faces[i] == -1:
            break
        # For each face, find the vertices of that face, and output it.
        # We make no effort to re-use vertices between multiple faces,
        # A fancier implementation might do so.
        edges = []
        for j in range(3):
            edges.append(faces[i:][j])
        verts = [edge_to_boundary_vertex(edge, f_eval, x, y, z) for edge in edges]
        next_vert_index = len(output_verts) + 1
        tri = Triangle(
            next_vert_index,
            next_vert_index+1,
            next_vert_index+2,
        )
        output_verts.extend(verts)
        output_tris.append(tri)
    return Mesh(output_verts, output_tris)


def marching_cubes_3d(f, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX, zmin=ZMIN, zmax=ZMAX):
    """Iterates over a cells of size one between the specified range, and evaluates f to produce
        a boundary by Marching Cubes. Returns a Mesh object."""
    # For each cube, evaluate independently.
    # If this wasn't demonstration code, you might actually evaluate them together for efficiency
    mesh = Mesh()
    for x in range(xmin, xmax):
        for y in range(ymin, ymax):
            for z in range(zmin, zmax):
                cell_mesh = marching_cubes_3d_single_cell(f, x, y, z)
                mesh.extend(cell_mesh)
    return mesh


def circle(x, y, z):
    return math.sqrt(x*x + y*y + z*z) - ISOVALUE

def octaeder(x, y, z):
    return abs(x) + abs(y) + abs(z) - ISOVALUE

def torus(x, y, z):
    return (x**2 + y**2 + z**2 +0.375)**2 - 2*(x**2 + y**2)

def cube(x, y, z):
    return max(abs(x), abs(y), abs(z)) - ISOVALUE

if __name__ == "__main__":
    mesh = marching_cubes_3d(cube)
    fileName = "cube"
    mesh.export_in_off(fileName)
    mesh.export_in_obj(fileName)
