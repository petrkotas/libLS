import numpy
import struct
import sys

octahedron_vertices = numpy.array( [ 
    [ 1.0, 0.0, 0.0], # 0 
    [-1.0, 0.0, 0.0], # 1
    [ 0.0, 1.0, 0.0], # 2 
    [ 0.0,-1.0, 0.0], # 3
    [ 0.0, 0.0, 1.0], # 4 
    [ 0.0, 0.0,-1.0]  # 5                                
    ] )
octahedron_triangles = numpy.array( [ 
    [ 0, 4, 2 ],
    [ 2, 4, 1 ],
    [ 1, 4, 3 ],
    [ 3, 4, 0 ],
    [ 0, 2, 5 ],
    [ 2, 1, 5 ],
    [ 1, 3, 5 ],
    [ 3, 0, 5 ]] )

ASCII_FACET = """facet normal 0 0 0
outer loop
vertex {face[0][0]:.4f} {face[0][1]:.4f} {face[0][2]:.4f}
vertex {face[1][0]:.4f} {face[1][1]:.4f} {face[1][2]:.4f}
vertex {face[2][0]:.4f} {face[2][1]:.4f} {face[2][2]:.4f}
endloop
endfacet
"""

BINARY_HEADER = "80sI"
BINARY_FACET = "12fH"

# class declarations

class Triangulated_Sphere:
    def normalize_v3(arr):
        ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
        lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
        arr[:,0] /= lens
        arr[:,1] /= lens
        arr[:,2] /= lens                
        return arr

    def divide_all( vertices, triangles ):    
        #new_triangles = []
        new_triangle_count = len( triangles ) * 4
        # Subdivide each triangle in the old approximation and normalize
        #  the new points thus generated to lie on the surface of the unit
        #  sphere.
        # Each input triangle with vertices labelled [0,1,2] as shown
        #  below will be turned into four new triangles:
        #
        #            Make new points
        #                 a = (0+2)/2
        #                 b = (0+1)/2
        #                 c = (1+2)/2
        #        1
        #       /\        Normalize a, b, c
        #      /  \
        #    b/____\ c    Construct new triangles
        #    /\    /\       t1 [0,b,a]
        #   /  \  /  \      t2 [b,1,c]
        #  /____\/____\     t3 [a,b,c]
        # 0      a     2    t4 [a,c,2]    
        v0 = vertices[ triangles[:,0] ]
        v1 = vertices[ triangles[:,1] ]
        v2 = vertices[ triangles[:,2] ]
        a = ( v0+v2 ) * 0.5
        b = ( v0+v1 ) * 0.5
        c = ( v1+v2 ) * 0.5  
        normalize_v3( a )
        normalize_v3( b )
        normalize_v3( c )
    
        #Stack the triangles together.
        vertices = numpy.vstack( (v0,b,a,  b,v1,c,  a,b,c, a,c,v2) )
        #Now our vertices are duplicated, and thus our triangle structure are unnecesarry.    
        return vertices, numpy.arange( len(vertices) ).reshape( (-1,3) )

    def create_unit_sphere( recursion_level=2 ):
        vertex_array, index_array = octahedron_vertices, octahedron_triangles
        for i in range( recursion_level - 1 ):
            vertex_array, index_array  = divide_all(vertex_array, index_array)
        return vertex_array, index_array


    def vertex_array_only_unit_sphere( recursion_level=2 ):
        vertex_array, index_array = create_unit_sphere(recursion_level)
        if recursion_level > 1:    
            return vertex_array.reshape( (-1) )
        else:
            return vertex_array[index_array].reshape( (-1) )





class ASCII_STL_Writer:
    """ Export 3D objects build of 3 or 4 vertices as ASCII STL file.
    """
    def __init__(self, stream):
        self.fp = stream
        self._write_header()

    def _write_header(self):
        self.fp.write("solid python\n")

    def close(self):
        self.fp.write("endsolid python\n")

    def _write(self, face):
        self.fp.write(ASCII_FACET.format(face=face))

    def _split(self, face):
        p1, p2, p3, p4 = face
        return (p1, p2, p3), (p3, p4, p1)

    def add_face(self, face):
        """ Add one face with 3 or 4 vertices. """
        if len(face) == 4:
            face1, face2 = self._split(face)
            self._write(face1)
            self._write(face2)
        elif len(face) == 3:
            self._write(face)
        else:
            raise ValueError('only 3 or 4 vertices for each face')

    def add_faces(self, faces):
        """ Add many faces. """
        for face in faces:
            self.add_face(face)

class Binary_STL_Writer(ASCII_STL_Writer):
    """ Export 3D objects build of 3 or 4 vertices as binary STL file.
    """
    def __init__(self, stream):
        self.counter = 0
        super(Binary_STL_Writer, self).__init__(stream)

    def close(self):
        self._write_header()

    def _write_header(self):
        self.fp.seek(0)
        self.fp.write(struct.pack(BINARY_HEADER, b'Python Binary STL Writer', self.counter))

    def _write(self, face):
        self.counter += 1
        data = [
            0., 0., 0.,
            face[0][0], face[0][1], face[0][2],
            face[1][0], face[1][1], face[1][2],
            face[2][0], face[2][1], face[2][2],
            0
        ]
        self.fp.write(struct.pack(BINARY_FACET, *data))


def example():
    def get_cube():
        # cube corner points
        s = 3.
        p1 = (0, 0, 0)
        p2 = (0, 0, s)
        p3 = (0, s, 0)
        p4 = (0, s, s)
        p5 = (s, 0, 0)
        p6 = (s, 0, s)
        p7 = (s, s, 0)
        p8 = (s, s, s)

        # define the 6 cube faces
        # faces just lists of 3 or 4 vertices
        return [
            [p1, p5, p7, p3],
            [p1, p5, p6, p2],
            [p5, p7, p8, p6],
            [p7, p8, p4, p3],
            [p1, p3, p4, p2],
            [p2, p6, p8, p4],
        ]

    with open('cube.stl', 'wb') as fp:
        writer = Binary_STL_Writer(fp)
        writer.add_faces(get_cube())
        writer.close()

if __name__ == '__main__':
    if len(sys.argv) == 2:
        rl = sys.argv[1]
    else:
        rl = 3

    tr = Triangulated_Sphere()

    vertex_array, index_array = tr.create_unit_sphere()

    print vertex_array

