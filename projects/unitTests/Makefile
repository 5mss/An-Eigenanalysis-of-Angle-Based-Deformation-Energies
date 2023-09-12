include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/catch2 -I../../ext/solvePoly
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES} -std=c++17
LDFLAGS    = ${LDFLAGS_COMMON} -std=c++17
MAIN			 = unitTests.cpp

SOURCES    = ../../ext/catch2/catch_amalgamated.cpp \
             ../../ext/solvePoly/poly34.cpp \
						 ../../src/util/COLLISION_UTIL.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/util/MIN_HEAP.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/PCG.cpp \
						 ../../src/util/DIAGONAL.cpp \
						 ../../src/util/STRAND_DIAGONAL.cpp \
						 ../../src/util/BLOCK_SPARSE_MATRIX3.cpp \
						 ../../src/util/BLOCK_DIAGONAL_MATRIX3.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/SPHERE.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_STRETCHING.cpp \
						 ../../src/Hyperelastic/Strand/SIN_BENDING.cpp \
						 ../../src/Hyperelastic/Strand/TAN_BENDING.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_F_BENDING.cpp \
						 ../../src/Hyperelastic/Shell/QUADRATIC_F_BENDING.cpp \
						 ../../src/Hyperelastic/Shell/ARAP.cpp \
						 ../../src/Hyperelastic/Shell/STRETCHING.cpp

include ../include_bottom.mk
