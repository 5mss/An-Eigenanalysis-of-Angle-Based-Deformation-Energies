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
#ifndef TRIANGLE_MESH_H
#define TRIANGLE_MESH_H
#define DISABLE_CCD

#ifndef DISABLE_CCD
#include "CCD_RESPONSE.h"
#include <SAFE_CCD.h>
#endif

#include "SETTINGS.h"
#include "Hyperelastic/Shell/STRETCHING.h"
#include "Hyperelastic/Shell/BENDING.h"
#include "Hyperelastic/Volume/VERTEX_FACE_COLLISION.h"
#include "Hyperelastic/Volume/EDGE_COLLISION.h"

#include <map>
#include <vector>
#include <set>


namespace HOBAK {

using namespace std;

class TRIANGLE_MESH
{
public:
  TRIANGLE_MESH() = default;
  TRIANGLE_MESH(const vector<VECTOR3>& restVertices, 
                const vector<VECTOR3I>& triangles);
  virtual ~TRIANGLE_MESH();

  // accessors
  const vector<VECTOR3>& vertices() const     { return _vertices; };
  vector<VECTOR3>& vertices()                 { return _vertices; };
  const vector<VECTOR3>& verticesOld() const  { return _verticesOld; };
  vector<VECTOR3>& verticesOld()              { return _verticesOld; };
  const vector<VECTOR3>& restVertices() const { return _restVertices; };
  vector<VECTOR3>& restVertices()             { return _restVertices; };
  const vector<VECTOR3I>& triangles() const   { return _triangles; };
  vector<VECTOR3I>& triangles()               { return _triangles; };
  const vector<VECTOR2I>& edges() const       { return _edges; };
  /*
  const vector<VECTOR4I>& vertexFaceCollisionTets() const    { return _vertexFaceCollisionTets; };
  vector<VECTOR4I>& vertexFaceCollisionTets()                { return _vertexFaceCollisionTets; };
  */
  const vector<REAL>& restOneRingAreas() const     { return _restOneRingAreas; };
  vector<REAL>& restOneRingAreas()                 { return _restOneRingAreas; };
  const VECTOR& restEdgeAreas() const        { return _restEdgeAreas; };
  const vector<REAL>& restThetas() const           { return _restThetas; };
  vector<REAL>& restThetas()                       { return _restThetas; };
  const vector<VECTOR4I>& flaps() const            { return _flaps; };
  vector<VECTOR4I>& flaps()                        { return _flaps; };
  const VECTOR3& vertex(const int index) const     { return _vertices[index]; };
  VECTOR3& vertex(const int index)                 { return _vertices[index]; };
  const VECTOR3& vertexOld(const int index) const  { return _verticesOld[index]; };
  VECTOR3& vertexOld(const int index)              { return _verticesOld[index]; };
  const VECTOR3& restVertex(const int index) const { return _restVertices[index]; };
  const REAL & collisionEps() const                { return _collisionEps; };
  const vector<bool>& invertedVertices() const     {return _invertedVertices; };
  const vector<VECTOR2I>& edgeTriangleNeighbors() const {return _edgeTriangleNeighbors; };
  /*
  const VECTOR4I& tet(const int tetIndex) const    { return _tets[tetIndex]; };
  const vector<int>& surfaceTets() const           { return _surfaceTets; };
  const vector<int>& surfaceVertices() const       { return _surfaceVertices; };
  const vector<VECTOR3I>& surfaceTriangles() const { return _surfaceTriangles; };
  const vector<VECTOR2I>& surfaceEdges() const     { return _surfaceEdges; };
  const vector<pair<int,int> >& vertexFaceCollisions() const { return _vertexFaceCollisions; };
  const vector<pair<int,int> >& edgeEdgeCollisions() const   { return _edgeEdgeCollisions; };
  const vector<bool>& edgeEdgeIntersections() const          { return _edgeEdgeIntersections; };
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinates() const   { return _edgeEdgeCoordinates; };
  const vector<REAL>& surfaceTriangleAreas() const { return _surfaceTriangleAreas; };
  const vector<VECTOR3I>& surfaceTriangleNeighbors() const { return _surfaceTriangleNeighbors; };
  */

  int totalVertices() const { return _vertices.size(); };
  const int DOFs() const    { return _vertices.size() * 3; };

  // get deformation gradient, and its SVD
  MATRIX3x2 computeF(const int index) const;
  void computeFs();
  //void computeFdots(const VECTOR& velocity);
  void computeSVDs();

  // get volume-weighted global translation
  VECTOR3 getTranslation() const;
  
  // get volume-weighted global translation, for the rest state
  VECTOR3 getRestTranslation() const;
 
 /* 
  // get Procrustes-style global rotation
  MATRIX3 getRotation() const;
  */

  // get the current displacement in vector form
  VECTOR getDisplacement() const;
  
  // set the vertex displacements to these values exactly
  void setDisplacement(const VECTOR& delta);

  // set the vertex positions directly exactly
  void setPositions(const VECTOR& positions);

  // add the following deltas to the positions
  void addDisplacement(const VECTOR& delta);
  
  // set collision eps to something new
  void setCollisionEps(const REAL& eps);
  void setCollisionStiffness(const REAL& stiffness);
 
  /*
  // set collision pairs (for replays)
  void setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                         const vector<pair<int,int> >& edgeEdge);
  */

  // compute stretching quantities
  REAL computeStretchingEnergy(const SHELL::STRETCHING& hyperelastic) const;
  VECTOR computeStretchingForces(const SHELL::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingClampedHessian(const SHELL::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingHessian(const SHELL::STRETCHING& hyperelastic) const;

  REAL computeBendingEnergy(const SHELL::BENDING& hyperelastic) const;
  VECTOR computeBendingForces(const SHELL::BENDING& hyperelastic) const;
  SPARSE_MATRIX computeBendingHessian(const SHELL::BENDING& hyperelastic) const;
  SPARSE_MATRIX computeBendingClampedHessian(const SHELL::BENDING& hyperelastic) const;

  /*
  // compute damping quantities
  VECTOR computeDampingForces(const VOLUME::DAMPING& damping) const;
  virtual SPARSE_MATRIX computeDampingHessian(const VOLUME::DAMPING& damping) const;

  // compute collision quantities
  //VECTOR computeCollisionForces(const REAL& collisionStiffness) const;
  //SPARSE_MATRIX computeCollisionClampedHessian(const REAL& collisionStiffness) const;
  */

  // compute x-based collision quantities
  VECTOR computeVertexFaceCollisionForces() const;
  
  
  SPARSE_MATRIX computeVertexFaceCollisionClampedHessian() const;
  //VECTOR computeEberleDampingForces(const VECTOR& velocity, const REAL& integratorConstant,
  //                                  const REAL& collisionStiffness) const;
  //SPARSE_MATRIX computeEberleClampedHessian(const VECTOR& velocity, const REAL& integratorConstant,
  //                                          const REAL& collisionStiffness) const;

  REAL computeEdgeEdgeCollisionEnergy() const;
  VECTOR computeEdgeEdgeCollisionForces() const;
  //VECTOR computeEdgeEdgeDampingForces(const VECTOR& velocity, const REAL& integratorConstant,
  //                                    const REAL& collisionStiffness) const;
  SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian() const;
  //SPARSE_MATRIX computeEdgeEdgeDampingClampedHessian(const VECTOR& velocity, const REAL& integratorConstant,
  //                                                   const REAL& collisionStiffness) const;
    
  // Does chronological resolution to CCDs detected in edge-edge and vertex-face
  // and gives updated displacement from rest pose
  VECTOR computeCCDResponsePosition(const int maxIters);

  // Takes all the ee and vf CCDs and gathers them into _allCCDs in chronological order
  // updates _allCCDs
  void gatherCCDs();

  // Takes all CCDs and groups them into contiguous regions
  // updates _regionsCCD
  void gatherCCDRegions();

  // Loops through each CCD region and resolves the most current CCD
  void solveCCDRegions();

  /*
  // compuate elastic and damping forces at the same time
  virtual VECTOR computeInternalForces(const VOLUME::HYPERELASTIC& hyperelastic,
                                       const VOLUME::DAMPING& damping) const;

  // get the bounding box for the current mesh
  void getBoundingBox(VECTOR3& mins, VECTOR3& maxs) const;  
  */
  
  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions();

  // find CCD'd vertex face pairs using _vertices and _verticesOld
  virtual void computeVertexFaceCCD(); 
  
  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions();

  // find CCD'd edge-edge pairs using _vertices and _verticesOld
  virtual void computeEdgeEdgeCCD();

  // debug edge-edge collisions, load up some specific pairs
  //void computeEdgeEdgeCollisionsDebug();

  // based on vertex-face collision pairs, build "collision tets"
  void buildVertexFaceCollisionTets(const VECTOR& velocity);
  
  /*  
  // write out the surface to OBJ triangle mesh
  static bool writeSurfaceToObj(const string& filename, const TRIANGLE_MESH& tetMesh);
  */

  // read in OBJ as a triangle mesh
  static bool readObjFile(const string& filename, 
                           vector<VECTOR3>& vertices, vector<VECTOR3I>& triangles);

  
  // write out OBJ mesh file, the bool tells whether to write the
  // rest or deformed vertices
  static bool writeObjFile(const string& filename, const TRIANGLE_MESH& triMesh, const bool restVertices);
  
  // normalize vertices so that they're in a unit box, centered at (0.5, 0.5, 0.5)
  static vector<VECTOR3> normalizeVertices(const vector<VECTOR3>& vertices);
  
  // compute distance between a point and triangle
  static REAL pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                    const VECTOR3& v2, const VECTOR3& v);
  
  // see if the projection of v onto the plane of v0,v1,v2 is inside the triangle
  // formed by v0,v1,v2
  static bool pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                          const VECTOR3& v2, const VECTOR3& v);

  // see if a current surface triangle has been crushed to degeneracy
  bool surfaceTriangleIsDegenerate(const int surfaceTriangleID);

  // see if the vertex is inside the collision cell described in 
  // Chapter 11: Collision Processing, [Kim and Eberle 2020]
  bool insideCollisionCell(const int TriangleID, const VECTOR3& vertex);
  /*
  // compute the dihedral angle between two surface faces
  REAL surfaceFaceDihedralAngle(const int surfaceID0, const int surfaceID1) const;
  */
  // get the normal to a plane, specified by three points
  static VECTOR3 planeNormal(const vector<VECTOR3>& plane);

protected:
  /*
  // compute volumes for each vertex one-ring -- works for rest
  // and deformed, just pass it _restVertices or _vertices
  vector<REAL> computeOneRingVolumes(const vector<VECTOR3>& vertices);
  */

  // compute material inverses for deformation gradient
  vector<MATRIX2> computeDmInvs();

  // derivatives with respect to shape functions
  static MATRIX3x2 computeDshape(const int i);

  // compute change-of-basis from deformation gradient F to positions, x
  vector<MATRIX6x9> computePFpxs();

  void computeEdges();
  void computeSurfaceAreas();
  void computeTriangleNeighbors();
  void computeEdgeTriangleNeighbors();
  void computeFlaps();

  /*
  // (DEACTIVATED) find out how close all the edges are initially
  void computeEdgeEdgeRestDistance();
  */

  // compute a triangle area
  static REAL triangleArea(const vector<VECTOR3>& triangle);


  // project point onto plane, specific by three points
  static VECTOR3 pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point);
  
  
  
  // compute distance to collision cell wall, where positive means inside
  // and negative means outside
  //REAL distanceToCollisionCellWall(const int triangleID, const VECTOR3& vertex);
  
  // compute whether one vertex is inside the vertex one right of another
  void computeSurfaceVertexOneRings();

  /*
  // build a consistent tet/flap ordering from two surface triangles
  VECTOR4I buildSurfaceFlap(const int surfaceID0, const int surfaceID1) const;

  // compute the normal of the surface triangle at _surfaceTriangles[triangleID];
  VECTOR3 surfaceTriangleNormal(const int triangleID) const;
  
  */
  
  
  // compute which vertices are attached to inverted tets
  void computeInvertedVertices();
  
  // are these two triangles neighbors?
  bool areTriangleNeighbors(const int id0, const int id1) const;

  // build a consistent tet/flap ordering from two triangles
  VECTOR4I buildFlap(const int triangle0, const int triangle1) const;

  // the core geometry
  vector<VECTOR3>    _vertices;
  vector<VECTOR3>    _verticesOld;
  vector<VECTOR3>    _restVertices;
  vector<VECTOR3I>   _triangles;

  vector<VECTOR4I>   _flaps;
  vector<REAL>       _restFlapAreas;
  vector<REAL>       _restThetas;

  /*
  // volumes, computed by computeTetVolumes and computeOneRingVolumes
  //
  // TODO: why aren't these just VECTORs?
  vector<REAL>  _restTetVolumes;
  vector<REAL>  _restOneRingVolumes;
  */
  vector<REAL>  _restOneRingAreas;
  VECTOR        _restEdgeAreas;

  // support for computing deformation gradient F
  vector<MATRIX2>    _DmInvs;

  // change-of-basis to go from deformation gradient (F) to positions (x)
  vector<MATRIX6x9> _pFpxs;

  // deformation gradients, and their SVDs
  vector<MATRIX3x2> _Fs;
  vector<MATRIX3x2> _Us;
  vector<VECTOR2> _Sigmas;
  vector<MATRIX2> _Vs;

  /*
  // velocity gradients
  vector<MATRIX3> _Fdots;

  // list of tets that are on the surface
  vector<int> _surfaceTets;

  // list of triangles that are on the surface
  // each triplet is ordered counter-clockwise, facing outwards
  // the VECTOR3I indexes into _vertices
  vector<VECTOR3I> _surfaceTriangles;
  */
  vector<REAL> _restTriangleAreas;
  vector<REAL> _triangleAreas;

  // for each triangle, what's the index of the neighboring triangles?
  vector<VECTOR3I> _triangleNeighbors;

  // list of edges on the surface
  // each pair is in sorted order, and index into _vertices
  vector<VECTOR2I> _edges;

  /*
  // list of vertices that are on the surface
  // indexes into _vertices
  vector<int> _surfaceVertices;
  */

  // for each _edge, what are the one or two neighboring triangles
  // in _triangles?
  vector<VECTOR2I> _edgeTriangleNeighbors;

  /*
  // for each pair of _surfaceEdges, what _collisionEps should we use? If they started
  // out closer than _collisionEps, then we need to set a smaller tolerance.
  //
  // the entries in the pair<int,int> are:
  //   unsigned int flat = edge[0] + edge[1] * _surfaceEdges.size();
  map<pair<unsigned int, unsigned int>, REAL> _edgeEdgeRestDistance;
  */
  // how close is considered to be in collision?
  REAL _collisionEps;
  
  // The ccd detector object
  // SAFE_CCD<REAL>* _ccd;

  // list of vertex-face collisions and CCDs
  // (second.)first indexes into _vertices
  // (second.)second indexes into _surfaceTriangles
  // first indexes into collision time for CCD (CCDs sorted by collision time)
  vector<pair<int, int> > _vertexFaceCollisions;
  vector<pair<REAL, pair<int, int> > > _vertexFaceCCDs;

  
  // list of edge-edge collision indices and CCDs
  // (second.)first indexes into _surfaceEdges
  // (second.)second indexes into _surfaceEdges
  // first indexes into collision time for CCD (CCDs sorted by collision time)
  vector<pair<int, int> > _edgeEdgeCollisions;
  vector<pair<REAL, pair<int,int> > > _edgeEdgeCCDs;

  // .first represents time of collision
  // .second represents the index into _vertices with 5th int as type (0 for vf, 1 for ee)
  // sorted chronologically at population time in gatherCCDs()
  vector<pair<REAL,VECTOR5I>> _allCCDs;
  
  // Each element in overlying vector represents a contiguous region. 
  // Individual vectors are collections of indices into _allCCDs and should be sorted chrono
  vector<vector<int>> _regionsCCD;
  // Direct helper with _regionsCCD: each vector element is a zone
  // each individual vector contains the set of all indexes into _vertices that make up that zone
  vector<set<int>> _vertexRegionsCCD;

  // interpolation coordinates for edge-edge collisions
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinates;
  
  // are the edge-edge collisions still separate, or is there already a face-edge intersection?
  vector<bool> _edgeEdgeIntersections;

  // list of collisionEps for each _edgeEdgeIntersections. Usually this will be the default,
  // but if it started out closer than that, we use that instead
  //vector<REAL> _edgeEdgeCollisionEps;

  // list of "collision tets" formed by vertex-face pairs
  vector<VECTOR4I> _vertexFaceCollisionTets;
  
  /*
  // list of "collision tets" formed by edge-edge pairs
  //vector<VECTOR4I> _edgeEdgeCollisionsTets;

  // DEBUG: see if the collision tet exists already
  //map<pair<int, int>, int> _vertexFaceCollisionTetsHash;

  */

  // scaling term for vertex-face collision forces
  vector<REAL> _vertexFaceCollisionAreas;
  
  
  // scaling term for edge-edge collision forces
  vector<REAL> _edgeEdgeCollisionAreas;

  /*
  // convert into into _vertices into index into _surfaceVertices
  map<int, int> _volumeToSurfaceID;

  // constitutive model for collisions
  VOLUME::HYPERELASTIC* _collisionMaterial;
  */

  // have your computed the SVDs since the last time you computed F?
  bool _svdsComputed;

  // see if two indices in _vertices (in sorted order)
  // are within the one ring of each other
  map<pair<int,int>, bool> _insideSurfaceVertexOneRing;

  VOLUME::VERTEX_FACE_COLLISION* _vertexFaceEnergy;

  VOLUME::EDGE_COLLISION* _edgeEdgeEnergy;

  // which vertices are inverted?
  vector<bool> _invertedVertices;
};

}

#endif
