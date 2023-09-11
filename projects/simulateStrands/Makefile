include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/tinyobjloader -I../../ext/ -I/usr/local/include
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES}
#CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES} -DNDEBUG
LDFLAGS    = ${LDFLAGS_COMMON}

MAIN = simulateStrands.cpp
SOURCES    = ../../ext/glvu/glvu.cpp \
             ../../ext/solvePoly/poly34.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/util/COLLISION_UTIL.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/DRAW_GL.cpp \
						 ../../src/util/PCG.cpp \
						 ../../src/util/STRAND_DIAGONAL.cpp \
						 ../../src/util/DIAGONAL.cpp \
						 ../../src/util/BLOCK_SPARSE_MATRIX3.cpp \
						 ../../src/util/BLOCK_DIAGONAL_MATRIX3.cpp \
						 ../../src/Timestepper/Strand/TIMESTEPPER.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/CAPSULE.cpp \
						 ../../src/Geometry/CYLINDER.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/SPHERE.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Geometry/TRIANGLE_MESH.cpp \
						 ../../src/Hyperelastic/Volume/HYPERELASTIC.cpp \
						 ../../src/Hyperelastic/Volume/EDGE_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/VERTEX_FACE_COLLISION.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_STRETCHING.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_F_BENDING.cpp \
						 ../../src/Hyperelastic/Strand/TAN_BENDING.cpp \
						 ../../src/Geometry/STRAND_MESH.cpp \
						 ../../src/Geometry/STRAND_MESH_FASTER.cpp \
						 ../../src/Geometry/STRAND_NET_MESH.cpp

include ../include_bottom.mk
