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
#include "util/DRAW_GL.h"
#include <random>

///////////////////////////////////////////////////////////////////////
// A bunch fo drawing routines for ANGLE objects, all in one place
//
// This is the only file with a GL dependency, so to remove it,
// just don't make any of the draw calls listed here.
///////////////////////////////////////////////////////////////////////

#include <glvu.h>
#if __APPLE__
#include <GL/glut.h>
#elif __linux__
#include <GL/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif

#include <iostream>

namespace HOBAK {

using namespace std;

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  //glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);
  glOrtho(0,0,0,0, -10, 10);

  // set the matric mode back to modelview
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // must set color before setting raster position, otherwise it won't take
  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);

  // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
  //glRasterPos3f(-halfZoom* 0.95, -halfZoom* 0.95, 0);
  //glRasterPos3f(-0.5, -0.5, 0);
  //glRasterPos3f(-0.5, -0.5, 0);
  glRasterPos3f(-0.95, -0.95, 0);

  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawPlaneConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = triangleMesh->vertices();
  const vector<PLANE_CONSTRAINT>& constraints = stepper->planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = constraints[x];
    const KINEMATIC_SHAPE* shape = constraint.shape;

    int index = constraint.vertexID;
    VECTOR3 vertex = vertices[index];
    const VECTOR3& localClosest = constraint.localClosestPoint;
    const VECTOR3& localNormal = constraint.localNormal;

    VECTOR3 closestPoint = shape->localVertexToWorld(localClosest); 

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // connect to closest point
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 10.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // draw the normal too
    VECTOR3 normal = shape->localNormalToWorld(localNormal);
    normal *= 0.1;
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
      glVertex3f(closestPoint[0] + normal[0], closestPoint[1] + normal[1], closestPoint[2] + normal[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawPlaneConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper)
{
  //const vector<VECTOR3>& vertices = strandMesh->vertices();
  const vector<VECTOR3>& vertices = strandMesh->verticesOld();
  const vector<PLANE_CONSTRAINT>& constraints = stepper->planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = constraints[x];
    const KINEMATIC_SHAPE* shape = constraint.shape;

    int index = constraint.vertexID;
    VECTOR3 vertex = vertices[index];
    const VECTOR3& localClosest = constraint.localClosestPoint;
    const VECTOR3& localNormal = constraint.localNormal;

    VECTOR3 closestPoint = shape->localVertexToWorld(localClosest); 

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // connect to closest point
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // draw the normal too
    VECTOR3 normal = shape->localNormalToWorld(localNormal);
    normal *= 0.1;
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
      glVertex3f(closestPoint[0] + normal[0], closestPoint[1] + normal[1], closestPoint[2] + normal[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawKinematicConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = triangleMesh->vertices();
  const vector<KINEMATIC_CONSTRAINT>& constraints = stepper->kinematicConstraints();

  glPointSize(10.0);

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = constraints[x];
    int index = constraint.vertexID;

    const VECTOR3& vertex = vertices[index];
    glBegin(GL_POINTS);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawKinematicConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = strandMesh->vertices();
  const vector<KINEMATIC_CONSTRAINT>& constraints = stepper->kinematicConstraints();

  glPointSize(10.0);

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = constraints[x];
    int index = constraint.vertexID;

    glDisable(GL_DEPTH_TEST);
    const VECTOR3& vertex = vertices[index];
    glBegin(GL_POINTS);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  }
}



///////////////////////////////////////////////////////////////////////
// just draw a single vertex of a triangle mesh
///////////////////////////////////////////////////////////////////////
void drawVertex(const TRIANGLE_MESH& mesh, const int index)
{
  const vector<VECTOR3>& vertices = mesh.vertices();

  if (index < 0) return;
  if (index >= (int)vertices.size()) return;

  glBegin(GL_POINTS);
    const VECTOR3& v = vertices[index];
    glVertex3dv(v.data());
  glEnd();
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
}


///////////////////////////////////////////////////////////////////////
// draw the triangles of a shell
///////////////////////////////////////////////////////////////////////
void drawTriangleMesh(const TRIANGLE_MESH& mesh, bool drawOutlines)
{
  const vector<VECTOR3I>& triangles = mesh.triangles();
  VECTOR3 v[3];

  const VECTOR3 triangleColor(0.5, 0.5, 0.5);
  const VECTOR3 outlineColor(0, 0, 0);

  // draw front-facing triangles
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glColor4f(triangleColor[0], triangleColor[1], triangleColor[2], 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();
 
  // draw back-facing a different color 
  glCullFace(GL_FRONT);
  glColor4f(1.0, 0.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();

  // see if we're done
  if (!drawOutlines) return;

  glCullFace(GL_BACK);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);
    glLineWidth(2.0);
    //glColor4f(0, 0, 0, 1.0);
    glColor4f(outlineColor[0], outlineColor[1], outlineColor[1], 1.0);
    glBegin(GL_LINE_STRIP);
      for (int y = 0; y < 4; y++)
        glVertex3f(v[y % 3][0], v[y % 3][1], v[y % 3][2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
// draw a single strand
///////////////////////////////////////////////////////////////////////
void drawStrand(const STRAND_MESH& mesh, const int highlight)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  if (highlight >= strandIndices.size())
    return;

  // draw the strands
  glColor4f(0,0,0, 1.0);
  glLineWidth(5.0);
  glBegin(GL_LINES);
    for (unsigned int x = 0; x < strandIndices[highlight].size() - 1; x++)
    {
      const int i0 = strandIndices[highlight][x];
      const int i1 = strandIndices[highlight][x + 1];
      const VECTOR3& v0 = vertices[i0];
      const VECTOR3& v1 = vertices[i1];
      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
    }
  glEnd();

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    glColor4f(1.0, 0.0, 0.0, 1.0);
    glNormal3f(1.0, 0,0);
    for (unsigned int x = 0; x < strandIndices[highlight].size(); x++)
    {
      const int i0 = strandIndices[highlight][x];
      const VECTOR3& v0 = vertices[i0];
      glVertex3f(v0[0], v0[1], v0[2]);
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawStrandMesh(const STRAND_MESH& mesh, const int highlight, const bool drawVertices)
{
  // cout<<"----------in"<<endl;
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  // draw the strands
  glColor4f(0,0,0, 1.0);
  // glLineWidth(5.0);
  glLineWidth(2.0);
  glBegin(GL_LINES);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
      for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
      {
        const int i0 = strandIndices[y][x];
        const int i1 = strandIndices[y][x + 1];
        const VECTOR3& v0 = vertices[i0];
        const VECTOR3& v1 = vertices[i1];
        glVertex3f(v0[0], v0[1], v0[2]);
        glVertex3f(v1[0], v1[1], v1[2]);
      }
  glEnd();

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}



///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawStrandMeshOld(const STRAND_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.verticesOld();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  std::mt19937 gen(123456);
  //std::uniform_real_distribution<double> twister(0.0, 1.0);
  std::uniform_real_distribution<double> twister(0.0, 0.25);

  // draw the strands
  glLineWidth(1.0);
  glBegin(GL_LINES);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      // color each strand randomly
      REAL r = twister(gen);
      REAL g = twister(gen);
      REAL b = twister(gen);
      //glColor4f(r,g,b, 0.5);
      glColor4f(r,r,r,1.0);
      for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
      {
        const int i0 = strandIndices[y][x];
        const int i1 = strandIndices[y][x + 1];
        const VECTOR3& v0 = vertices[i0];
        const VECTOR3& v1 = vertices[i1];
        glVertex3f(v0[0], v0[1], v0[2]);
        glVertex3f(v1[0], v1[1], v1[2]);
      }
    }
  glEnd();

  if (!drawVertices) return;
  
  //glPointSize(10.0);
  glPointSize(3.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}



///////////////////////////////////////////////////////////////////////
// draw a sphere
///////////////////////////////////////////////////////////////////////
void drawSphere(const SPHERE& sphere)
{
  const MATRIX3& S = sphere.scale();
  const VECTOR3& t = sphere.translation();
  const Eigen::AngleAxis<GLfloat> R{ sphere.rotation().cast<GLfloat>() };

  //glColor4f(1.0, 0.0, 0.0, 0.5);
  glColor4f(0.99, 0.99, 0.99, 0.9);
  glPushMatrix();
    glTranslatef(t[0], t[1], t[2]);
    glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());
    glScalef(S(0,0), S(1,1), S(2,2));
    //glutSolidSphere(1.0, 20, 20);
    glutSolidSphere(1.0, 100, 100);
  glPopMatrix();
}


///////////////////////////////////////////////////////////////////////
// draw an AABB
///////////////////////////////////////////////////////////////////////
void drawAABB(const VECTOR3& minCorner, const VECTOR3& maxCorner)
{
  const VECTOR3 v000(minCorner[0], minCorner[1], minCorner[2]); 
  const VECTOR3 v100(maxCorner[0], minCorner[1], minCorner[2]); 
  const VECTOR3 v010(minCorner[0], maxCorner[1], minCorner[2]); 
  const VECTOR3 v110(maxCorner[0], maxCorner[1], minCorner[2]); 
  const VECTOR3 v001(minCorner[0], minCorner[1], maxCorner[2]); 
  const VECTOR3 v101(maxCorner[0], minCorner[1], maxCorner[2]); 
  const VECTOR3 v011(minCorner[0], maxCorner[1], maxCorner[2]); 
  const VECTOR3 v111(maxCorner[0], maxCorner[1], maxCorner[2]); 

  glColor4f(0.0, 1.0, 0.0, 0.25);
  glBegin(GL_QUADS);
    // x plus
    VECTOR3 normal = (v000 - v100).cross(v000 - v110);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v010.data());
    glVertex3dv(v110.data());
    glVertex3dv(v100.data());
    glVertex3dv(v000.data());

    // x minus
    normal = (v001 - v101).cross(v001 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v101.data());
    glVertex3dv(v111.data());
    glVertex3dv(v011.data());

    // y minus
    normal = (v000 - v100).cross(v000 - v101);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v000.data());
    glVertex3dv(v100.data());
    glVertex3dv(v101.data());
    glVertex3dv(v001.data());

    // y plus
    normal = (v010 - v110).cross(v010 - v111);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v011.data());
    glVertex3dv(v111.data());
    glVertex3dv(v110.data());
    glVertex3dv(v010.data());

    // z plus
    normal = (v000 - v010).cross(v000 - v011);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v011.data());
    glVertex3dv(v010.data());
    glVertex3dv(v000.data());

    // z minus
    normal = (v100 - v110).cross(v100 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v100.data());
    glVertex3dv(v110.data());
    glVertex3dv(v111.data());
    glVertex3dv(v101.data());
  glEnd();
}
void drawAABB(const AABB_NODE& node)
{
  drawAABB(node.mins, node.maxs);
}

///////////////////////////////////////////////////////////////////////
// draw a AABB tree at a specific depth
///////////////////////////////////////////////////////////////////////
void drawAABBTree(const AABB_NODE* node, const int drawDepth, const int currentDepth)
{
  if (node == NULL) return;

  if (currentDepth == drawDepth)
  {
    drawAABB(*node);
    return;
  }

  drawAABBTree(node->child[0], drawDepth, currentDepth + 1);
  drawAABBTree(node->child[1], drawDepth, currentDepth + 1);
}

///////////////////////////////////////////////////////////////////////
// draw a AABB tree at a specific depth
///////////////////////////////////////////////////////////////////////
void drawAABBTree(const AABB_TREE& tree, const int drawDepth)
{
  drawAABBTree(&(tree.root()), drawDepth, 0);
}

///////////////////////////////////////////////////////////////////////
// draw a cube
///////////////////////////////////////////////////////////////////////
void drawCube(const CUBE& cube)
{
  const MATRIX3& S = cube.scale();
  const MATRIX3& R = cube.rotation();
  const VECTOR3& t = cube.translation();
  const VECTOR3& minCorner = S * VECTOR3(-0.5, -0.5, -0.5);
  const VECTOR3& maxCorner = S * VECTOR3(0.5, 0.5, 0.5);

  VECTOR3 v000(minCorner[0], minCorner[1], minCorner[2]); 
  VECTOR3 v100(maxCorner[0], minCorner[1], minCorner[2]); 
  VECTOR3 v010(minCorner[0], maxCorner[1], minCorner[2]); 
  VECTOR3 v110(maxCorner[0], maxCorner[1], minCorner[2]); 
  VECTOR3 v001(minCorner[0], minCorner[1], maxCorner[2]); 
  VECTOR3 v101(maxCorner[0], minCorner[1], maxCorner[2]); 
  VECTOR3 v011(minCorner[0], maxCorner[1], maxCorner[2]); 
  VECTOR3 v111(maxCorner[0], maxCorner[1], maxCorner[2]); 

  v000 = R * v000 + t;
  v100 = R * v100 + t;
  v010 = R * v010 + t;
  v110 = R * v110 + t;
  v001 = R * v001 + t;
  v101 = R * v101 + t;
  v011 = R * v011 + t;
  v111 = R * v111 + t;

  glColor4f(0.8, 0.5, 0.0, 0.5);
  glBegin(GL_QUADS);
    // x plus
    VECTOR3 normal = (v000 - v100).cross(v000 - v110);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v010.data());
    glVertex3dv(v110.data());
    glVertex3dv(v100.data());
    glVertex3dv(v000.data());

    // x minus
    normal = (v001 - v101).cross(v001 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v101.data());
    glVertex3dv(v111.data());
    glVertex3dv(v011.data());

    // y minus
    normal = (v000 - v100).cross(v000 - v101);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v000.data());
    glVertex3dv(v100.data());
    glVertex3dv(v101.data());
    glVertex3dv(v001.data());

    // y plus
    normal = (v010 - v110).cross(v010 - v111);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v011.data());
    glVertex3dv(v111.data());
    glVertex3dv(v110.data());
    glVertex3dv(v010.data());

    // z plus
    normal = (v000 - v010).cross(v000 - v011);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v011.data());
    glVertex3dv(v010.data());
    glVertex3dv(v000.data());

    // z minus
    normal = (v100 - v110).cross(v100 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v100.data());
    glVertex3dv(v110.data());
    glVertex3dv(v111.data());
    glVertex3dv(v101.data());
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw coordinate axes, xyz = rgb
///////////////////////////////////////////////////////////////////////
void drawAxes()
{
  // draw coordinate axes
  glPushMatrix();
  //glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  // x axis is red
  glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
  glVertex3f(10.0f, 0.0f, 0.0f);

  // y axis is green
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 10.0f, 0.0f);

  // z axis is blue
  glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 10.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}


///////////////////////////////////////////////////////////////////////
void drawKinematicShape(const KINEMATIC_SHAPE& shape)
{
  using namespace std;

  const string& name = shape.name();

  if (name.compare(string("CUBE")) == 0)
  {
    drawCube((const CUBE&)shape);
    return;
  }

  if (name.compare(string("SPHERE")) == 0)
  {
    drawSphere((const SPHERE&)shape);
    return;
  }
  
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawCollisions(const STRAND_MESH& mesh)
{
  const vector<pair<int, int> >& collisions = mesh.edgeEdgeCollisions();
  const vector<VECTOR2I>& edges = mesh.edgeIndices();
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<pair<VECTOR2, VECTOR2> >& coordinates = mesh.edgeEdgeCoordinates();

  glDisable(GL_DEPTH_TEST);

  glLineWidth(10.0);
  glColor4f(0.0, 0.0, 1.0, 100.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glColor4f(1.0, 1.0, 0.0, 100.0);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;
    
    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glEnable(GL_DEPTH_TEST);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawCollisionsOld(const STRAND_MESH& mesh)
{
  const vector<pair<int, int> >& collisions = mesh.edgeEdgeCollisionsOld();
  const vector<VECTOR2I>& edges = mesh.edgeIndices();
  const vector<VECTOR3>& vertices = mesh.verticesOld();
  const vector<pair<VECTOR2, VECTOR2> >& coordinates = mesh.edgeEdgeCoordinatesOld();

  glDisable(GL_DEPTH_TEST);

  glLineWidth(10.0);
  glColor4f(0.0, 0.0, 1.0, 100.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glColor4f(1.0, 1.0, 0.0, 100.0);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    //cout << " coordinate 0: " << coordinate0.transpose() << endl;
    //cout << " coordinate 1: " << coordinate1.transpose() << endl;
    
    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];

    //const VECTOR3 diff = middle0 - middle1;
    //cout << " collision diff: " << diff.norm() << endl;
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glEnable(GL_DEPTH_TEST);
}


} // ANGLE
