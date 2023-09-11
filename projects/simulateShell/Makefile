include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/tinyobjloader
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES}
LDFLAGS    = ${LDFLAGS_COMMON}

MAIN = simulateShell.cpp
SOURCES    = ../../ext/glvu/glvu.cpp \
             ../../ext/solvePoly/poly34.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/util/COLLISION_UTIL.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/DRAW_GL.cpp \
						 ../../src/util/BLOCK_DIAGONAL_MATRIX3.cpp \
						 ../../src/util/BLOCK_SPARSE_MATRIX3.cpp \
						 ../../src/Timestepper/Shell/TIMESTEPPER.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/CAPSULE.cpp \
						 ../../src/Geometry/CYLINDER.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/SPHERE.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Geometry/TRIANGLE_MESH.cpp \
						 ../../src/Geometry/TRIANGLE_MESH_FASTER.cpp \
						 ../../src/Hyperelastic/Volume/HYPERELASTIC.cpp \
						 ../../src/Hyperelastic/Volume/EDGE_COLLISION.cpp \
						 ../../src/Hyperelastic/Volume/VERTEX_FACE_COLLISION.cpp \
						 ../../src/Hyperelastic/Shell/QUADRATIC_F_BENDING.cpp \
						 ../../src/Hyperelastic/Shell/ARAP.cpp \
						 ../../src/Hyperelastic/Shell/STRETCHING.cpp

include ../include_bottom.mk
