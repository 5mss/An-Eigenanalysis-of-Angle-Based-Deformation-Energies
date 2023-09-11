/*
This file is part of HOBAK.

HOBAK is free software: you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

HOBAK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with HOBAK. 
If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef DRAW_GL_H
#define DRAW_GL_H

#include "Geometry/CAPSULE.h"
#include "Geometry/CYLINDER.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/AABB_TREE.h"
#include "util/COLLISION_UTIL.h"
#include <glvu.h>

// DEBUG: here while TRIANGLE_MESH and SHELL::TIMESTEPPER are prototyping classes
#include "Geometry/STRAND_MESH.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"

namespace HOBAK {

// Print a string to the GL window
void printGlString(string output);

VECTOR3 planeNormal(const vector<VECTOR3>& plane);

// DEBUG: here while TRIANGLE_MESH and SHELL::TIMESTEPPER are prototyping classes
void drawKinematicConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper);
void drawPlaneConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper);
void drawKinematicConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper);
void drawPlaneConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper);

// draw the collision tets that have been built for this mesh

// just draw a single vertex of a mesh
void drawVertex(const TRIANGLE_MESH& mesh, const int index);

// draw only the surface triangles of a tet mesh
void drawTriangleMesh(const TRIANGLE_MESH& mesh, bool drawOutlines);

// draw strands, or strand information
void drawStrandMesh(const STRAND_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawStrandMeshOld(const STRAND_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawStrand(const STRAND_MESH& mesh, const int highlight = -1);
void drawCollisions(const STRAND_MESH& mesh);
void drawCollisionsOld(const STRAND_MESH& mesh);

// draw lines between differences between the two meshes

// draw kinematic objects
void drawKinematicShape(const KINEMATIC_SHAPE& shape);
void drawCapsule(const CAPSULE& capsule);
void drawCylinder(const CYLINDER& cylinder);
void drawSphere(const SPHERE& sphere);
void drawCube(const CUBE& cube);

// draw an AABB
void drawAABB(const VECTOR3& minCorner, const VECTOR3& maxCorner);
void drawAABB(const AABB_NODE& node);

// draw a AABB tree at a specific depth
void drawAABBTree(const AABB_NODE* node, const int drawDepth, const int currentDepth);
void drawAABBTree(const AABB_TREE& tree, const int drawDepth);

// draw coordinate axes, xyz = rgb
void drawAxes();


} // HOBAK

#endif
