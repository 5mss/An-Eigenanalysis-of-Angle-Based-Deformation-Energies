#include <cmath>
#include <iostream>
#include <float.h>

#define CATCH_CONFIG_MAIN

#include "util/MATRIX_UTIL.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/QUADRATIC_F_BENDING.h"
#include "Hyperelastic/Strand/TAN_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_F_BENDING.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/STRAND_MESH_FASTER.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/LINE_INTERSECT.h"
#include "Timestepper/Strand/TIMESTEPPER.h"
#include "util/COLLISION_UTIL.h"
#include "catch_amalgamated.hpp"
#include "ext/solvePoly/poly34.h"

using namespace HOBAK;
using namespace std;

#include "STRAND_TESTS.h"
#include "SHELL_TESTS.h"

