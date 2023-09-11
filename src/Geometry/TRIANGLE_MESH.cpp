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
#include "TRIANGLE_MESH.h"
#include "util/MATRIX_UTIL.h"
#include "util/COLLISION_UTIL.h"
#include "LINE_INTERSECT.h"
#include "util/TIMER.h"
#include <iostream>
#include <float.h> // need for FLT_MAX

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#define VERY_VERBOSE 0

// Alvin's project: toggle usage of collision energy: 
// 0 for geometry-based sqrt(vf)
// 1 for F-based sqrt(vf)
// 2 for mcAdam's
// 3 for vf
// 4 for F-based vf
// 5 for Adaptive spring constant (F-based wrong for now)
#define FABRIC_USE_COLLISION 0

// Alvin's project: toggling usage of edge collision energies:
// 0 for geometry-based edge-edge
// 1 for F-based edge-edge
// 2 for F-based sqrt(ee)
// 3 for geometry-based sqrtee
// 4 for geometry-bsaed hybrid
#define FABRIC_USE_EDGE 0

namespace HOBAK {

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const vector<VECTOR3>& restVertices, 
                             const vector<VECTOR3I>& triangles) :
    _vertices(restVertices),
    _verticesOld(restVertices),
    _restVertices(restVertices),
    _triangles(triangles)
{
  /*
  _restTetVolumes = computeTetVolumes(_restVertices);
  _restOneRingVolumes = computeOneRingVolumes(_restVertices);
  */
  _DmInvs = computeDmInvs();
  _pFpxs = computePFpxs();

  const int totalTriangles = _triangles.size();
  _Fs.resize(totalTriangles);
  _Us.resize(totalTriangles);
  _Sigmas.resize(totalTriangles);
  _Vs.resize(totalTriangles);

  //computeSurfaceTriangles();
  //computeSurfaceVertices();
  computeEdges();
  computeSurfaceAreas();
  computeTriangleNeighbors();
  computeEdgeTriangleNeighbors();

  computeFlaps();

  
  // one centimeter
  // _collisionEps = 0.01;
  // two centimeters, one seems to get into trouble without CCD
  _collisionEps = 0.02;
  
  //_collisionMaterial = NULL;  // this is still experimental
  _svdsComputed = false;

  // this gets overwritten by TIMESTEPPER every step, so a dummy is fine
  REAL stiffness = 1000.0;
  
  // must be called after _collisionEps has been set
  //computeEdgeEdgeRestDistance();
  
  // store which surface vertices are within the one rings of
  // each other
  computeSurfaceVertexOneRings();

  // if you want to try out the McAdams energy, 
  // here's the place to swap it in
  //_vertexFaceEnergy = new VOLUME::VERTEX_FACE_COLLISION(stiffness, _collisionEps); // default
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
 
  _vertexFaceEnergy = new VOLUME::VERTEX_FACE_COLLISION(stiffness, _collisionEps);
  cout << "Triangle mesh using " << _vertexFaceEnergy->name() << endl;
  //_vertexFaceEnergy = new VOLUME::MU_LINEAR_FILTER(stiffness, _collisionEps);
  //_vertexFaceEnergy = new VOLUME::MU_CUBIC_FILTER(stiffness, _collisionEps);

  //FOR NOW: See top of file to change which edge energy
  #if FABRIC_USE_EDGE == 0
  _edgeEdgeEnergy = new VOLUME::EDGE_COLLISION(stiffness, _collisionEps); // default
  #elif FABRIC_USE_EDGE == 1
  _edgeEdgeEnergy = new VOLUME::EDGE_COLLISION_F(stiffness, _collisionEps);
  #elif FABRIC_USE_EDGE == 2
  _edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION_F(stiffness, _collisionEps);
  #elif FABRIC_USE_EDGE == 3
  _edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION(stiffness, _collisionEps);
  #elif FABRIC_USE_EDGE == 4
  _edgeEdgeEnergy = new VOLUME::EDGE_HYBRID_COLLISION(stiffness, _collisionEps);
  #endif
  cout << "Triangle mesh edge using " << _edgeEdgeEnergy->name() << endl;
  
  // _ccd = new SAFE_CCD<REAL>();
  // _ccd->Set_Coefficients(0.25);

  // verified on the bunny drop scene 
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
  //_edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION(stiffness, _collisionEps);
  
  // verified on the bunny drop scene 
  //_vertexFaceEnergy = new VOLUME::MCADAMS_COLLISION(stiffness, _collisionEps);
  //_edgeEdgeEnergy = new VOLUME::EDGE_HYBRID_COLLISION(stiffness, _collisionEps);
  
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::~TRIANGLE_MESH()
{
  delete _vertexFaceEnergy;
  delete _edgeEdgeEnergy;
}

#if 0

///////////////////////////////////////////////////////////////////////
// compute rest volumes for each vertex one-ring
///////////////////////////////////////////////////////////////////////
vector<REAL> TRIANGLE_MESH::computeOneRingVolumes(const vector<VECTOR3>& vertices)
{
  const vector<REAL> tetVolumes = computeTetVolumes(vertices);
  unsigned int size = vertices.size();

  vector<REAL> oneRingVolumes(size);
  for (unsigned int x = 0; x < size; x++)
    oneRingVolumes[x] = 0.0;

  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    const REAL quarter = 0.25 * tetVolumes[x];
    for (int y = 0; y < 4; y++)
      oneRingVolumes[_tets[x][y]] += quarter;
  }

  return oneRingVolumes;
}
#endif

///////////////////////////////////////////////////////////////////////
// compute material inverses for deformation gradient
///////////////////////////////////////////////////////////////////////
vector<MATRIX2> TRIANGLE_MESH::computeDmInvs()
{
  vector<MATRIX2> DmInvs(_triangles.size());

  for (unsigned int f = 0; f < _triangles.size(); f++)
  {
    const VECTOR3I& face = _triangles[f];
    const VECTOR3& v0 = _restVertices[face[0]];
    const VECTOR3& v1 = _restVertices[face[1]];
    const VECTOR3& v2 = _restVertices[face[2]];

    MATRIX3x2 Dm;
    Dm.col(0) = v1 - v0;
    Dm.col(1) = v2 - v0;

    // compute the QR factorization
    MATRIX3x2 Qm;
    Qm.col(0) = Dm.col(0).normalized();
    Qm.col(1) = (Dm.col(1) - Dm.col(1).dot(Qm.col(0)) * Qm.col(0)).normalized();

    MATRIX2 Rm;
    Rm = Qm.transpose() * Dm;

    MATRIX2 RmInv;
    const double eps = 1e-8;

    if (fabs(Rm.determinant()) > eps)
      RmInv = Rm.inverse();
    else
    {
      Eigen::JacobiSVD<MATRIX2> svd(Rm, Eigen::ComputeThinU | Eigen::ComputeThinV);
      VECTOR2 singularInvs;
      singularInvs.setZero();
      for (unsigned int x = 0; x < 2; x++)
        if (svd.singularValues()[x] > 1e-8)
          singularInvs[x] = 1.0 / svd.singularValues()[x];

      RmInv = svd.matrixV() * singularInvs.asDiagonal() * svd.matrixU().transpose();
    }
    MATRIX2 DmInv = RmInv;
    DmInvs[f] = DmInv;
  }
  assert(DmInvs.size() == _triangles.size());
  return DmInvs;
}


//////////////////////////////////////////////////////////////////////////////
// derivatives with respect to shape functions
//////////////////////////////////////////////////////////////////////////////
MATRIX3x2 TRIANGLE_MESH::computeDshape(const int i)
{
  MATRIX3x2 dShape;
  dShape.setZero();

  if (i < 3)
  {
    dShape.row(i)[0] = -1;
    dShape.row(i)[1] = -1;
  }
  else if (i < 6)
    dShape.row(i - 3)[0] = 1;
  else
    dShape.row(i - 6)[1] = 1;

  return dShape;
}

///////////////////////////////////////////////////////////////////////
// compute change-of-basis from deformation gradient F to positions, x
///////////////////////////////////////////////////////////////////////
vector<MATRIX6x9> TRIANGLE_MESH::computePFpxs()
{
  vector<MATRIX6x9> pFpxs(_triangles.size());
  pFpxs.clear();
  for (unsigned int f = 0; f < _triangles.size(); f++)
  {
    MATRIX6x9 pFpu;

    MATRIX2 DmInv = _DmInvs[f];
    for (int i = 0; i < 9; i++)
    {
      MATRIX3x2 pDspu = computeDshape(i);
      MATRIX3x2 pFpuColumn = pDspu * DmInv;
      pFpu.col(i) = flatten(pFpuColumn);
    }
    pFpxs.push_back(pFpu);
  }

  return pFpxs;
}

#if 0
///////////////////////////////////////////////////////////////////////
// used by computeSurfaceTriangles as a comparator between two
// triangles
///////////////////////////////////////////////////////////////////////
struct triangleCompare
{
  bool operator()(const VECTOR3I &a, const VECTOR3I &b) const
  {
    if (a[0] < b[0])
      return true;
    if (a[0] > b[0])
      return false;

    if (a[1] < b[1])
      return true;
    if (a[1] > b[1])
      return false;

    if (a[2] < b[2])
      return true;
    if (a[2] > b[2])
      return false;
    return false;
  }
};

///////////////////////////////////////////////////////////////////////
// find which triangles are on the surface
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceTriangles()
{
  map<VECTOR3I, int, triangleCompare> faceCounts;

  // for each tet, add its faces to the face count
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    VECTOR4I t = _tets[x];

    VECTOR3I faces[4];
    faces[0] << t[0], t[1], t[3];
    faces[1] << t[0], t[2], t[1];
    faces[2] << t[0], t[3], t[2];
    faces[3] << t[1], t[2], t[3];

    for (int y = 0; y < 4; y++)
      std::sort(faces[y].data(), faces[y].data() + faces[y].size());

    for (int y = 0; y < 4; y++)
      faceCounts[faces[y]]++;
  }

  // go back through the tets, if any of its faces have a count less than
  // 2, then it must be because it faces outside
  _surfaceTriangles.clear();
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    VECTOR4I t = _tets[x];

    VECTOR3I faces[4];

    // these are consistently ordered counter-clockwise
    faces[0] << t[0], t[1], t[3];
    faces[1] << t[0], t[2], t[1];
    faces[2] << t[0], t[3], t[2];
    faces[3] << t[1], t[2], t[3];

    VECTOR3I facesSorted[4];

    // make a sorted copy, but keep the original around for rendering
    for (int y = 0; y < 4; y++)
    {
      facesSorted[y] = faces[y];
      std::sort(facesSorted[y].data(), facesSorted[y].data() + facesSorted[y].size());
    }

    // see which faces don't have a dual
    for (int y = 0; y < 4; y++)
    {
      if (faceCounts[facesSorted[y]] < 2)
        _surfaceTriangles.push_back(faces[y]);
    }
  }

  cout << " Found " << _surfaceTriangles.size() << " surface triangles out of " 
       << _tets.size() * 4 << " possible. " << endl;
}
#endif

///////////////////////////////////////////////////////////////////////
// for each edge, what're the indices of the neighboring triangles?
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeTriangleNeighbors()
{
  // translate the VEC2I into a index
  map<pair<int,int>, int> edgeToIndex;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> toHash;
    toHash.first  = _edges[x][0];
    toHash.second = _edges[x][1];

    if (toHash.first > toHash.second)
    {
      int temp = toHash.first;
      toHash.first = toHash.second;
      toHash.second = temp;
    }
    edgeToIndex[toHash] = x;
  }

  // look up the edges of each surface face, tabulate the adjacent
  // triangles
  vector<vector<int> > faceHash(_edges.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // store the edge
      assert(edgeToIndex.find(edge) != edgeToIndex.end());
      int edgeIndex = edgeToIndex[edge];
      faceHash[edgeIndex].push_back(i);
    }
  }

  // store the final results
  _edgeTriangleNeighbors.resize(_edges.size());
  for (unsigned int i = 0; i < _edges.size(); i++)
  {
    _edgeTriangleNeighbors[i][0] = -1;
    _edgeTriangleNeighbors[i][1] = -1;

    assert(faceHash[i].size() > 0);
    _edgeTriangleNeighbors[i][0] = faceHash[i][0];

    if (faceHash[i].size() == 2)
      _edgeTriangleNeighbors[i][1] = faceHash[i][1];
  }
}

///////////////////////////////////////////////////////////////////////
// for each triangle, what's the index of the neighboring triangles?
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeTriangleNeighbors()
{
  multimap<pair<int, int>, unsigned int> edgeNeighboringTriangles;

  // hash all the edges from each surface triangle
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // hash it
      pair<pair<int, int>, unsigned int> hash(edge, i);
      edgeNeighboringTriangles.insert(hash);
    }
  }
 
  // get the other edge that wasn't the current one 
  _triangleNeighbors.clear();
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store results here
    VECTOR3I neighbors(-1,-1,-1);

    // reconstruct the edge again
    pair<int,int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];
      // for cloth, we can have naked edges
      bool orphaned = true; 

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // find the matching triangles
      auto range = edgeNeighboringTriangles.equal_range(edge);
      for (auto it = range.first; it != range.second; it++)
      {
        if (it->second != i){
          neighbors[j] = it->second;
          orphaned = false;
        }
      }
      // if the edgeneighboringtriangles only finds the self triangle
      // input -(j+1).
      if (orphaned) neighbors[j] = -(j+1);
    }

    // store the neighbors
    _triangleNeighbors.push_back(neighbors);
  }
}

///////////////////////////////////////////////////////////////////////
// compute surface areas for collision weights
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceAreas()
{
  // compute the areas
  _triangleAreas.clear();
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    vector<VECTOR3> vertices(3);
    vertices[0] = _restVertices[_triangles[x][0]];
    vertices[1] = _restVertices[_triangles[x][1]];
    vertices[2] = _restVertices[_triangles[x][2]];

    _triangleAreas.push_back(triangleArea(vertices));
  }

  // cache these out for when triangles deform later
  _restTriangleAreas = _triangleAreas;
  
  // compute the one-ring areas
  assert(_vertices.size() != 0);
  _restOneRingAreas.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _restOneRingAreas[x] = 0;
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    assert(x < _triangleAreas.size());
    assert(x < _triangles.size());
    const REAL& area = _triangleAreas[x];
    const VECTOR3I& triangle = _triangles[x];

    for (int y = 0; y < 3; y++)
    {
      const int index = triangle[y];
      assert(index < (int)_restOneRingAreas.size());
      _restOneRingAreas[index] += (1.0 / 3.0) * area;
    }
  }

  // build a mapping from edge index pairs to _edges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    edgeHash[edge] = x;
  }

  // compute the edge areas
  assert(_edges.size() != 0);
  _restEdgeAreas.resize(_edges.size());
  _restEdgeAreas.setZero();
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // build each edge
    for (int y = 0; y < 3; y++)
    {
      pair<int,int> edge(_triangles[x][y],
                         _triangles[x][(y + 1) % 3]);

      // swap them to the order the hash expects
      if (edge.first > edge.second)
      {
        const REAL temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      const int edgeIndex = edgeHash[edge];
      assert(edgeIndex >= 0);
      assert(edgeIndex < _restEdgeAreas.size());
      _restEdgeAreas[edgeIndex] += _triangleAreas[x] / 3.0;
    }
  }
}

#if 0
///////////////////////////////////////////////////////////////////////
// find which vertices are on the surface
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceVertices()
{
  if (_surfaceTriangles.size() == 0)
    computeSurfaceTriangles();

  // hash them all out
  std::map<int, bool> foundVertices;
  for (unsigned int x = 0; x < _surfaceTriangles.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++)
      foundVertices[_surfaceTriangles[x][y]] = true;
  }

  // serialize
  _surfaceVertices.clear();
  for (auto iter = foundVertices.begin(); iter != foundVertices.end(); iter++)
    _surfaceVertices.push_back(iter->first);

  // compute the reverse lookup
  for (unsigned int x = 0; x < _surfaceVertices.size(); x++)
    _volumeToSurfaceID[_surfaceVertices[x]] = x;

  cout << " Found " << _surfaceVertices.size() << " vertices on the surface " << endl;
}
#endif

///////////////////////////////////////////////////////////////////////
// find all edges
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdges()
{
  // hash all the edges, so we don't store any repeats
  map<pair<int, int>, bool> edgeHash;
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++)
    {
      const int v0 = _triangles[x][y];
      const int v1 = _triangles[x][(y + 1) % 3];

      // store them in sorted order
      pair<int, int> edge;
      if (v0 > v1)
      {
        edge.first = v1;
        edge.second = v0;
      }
      else
      {
        edge.first = v0;
        edge.second = v1;
      }

      // hash it out
      edgeHash[edge] = true;
    }
  }

  // store all the unique hashes
  _edges.clear();
  for (auto iter = edgeHash.begin(); iter != edgeHash.end(); iter++)
  {
    const pair<int,int> e = iter->first;
    const VECTOR2I edge(e.first, e.second);
    _edges.push_back(edge);
  }

  cout << " Found " << _edges.size() << " edges on the surface " << endl;
}

///////////////////////////////////////////////////////////////////////
// get all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFs()
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Fs.size() == _triangles.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    _Fs[x] = computeF(x);
    /*
    if (x < 10)
    {
      cout << " F " << x  << ": " << endl;
      cout << _Fs[x] << endl;
    }
    */

    assert(!_Fs[x].hasNaN());
  }

  // SVDs are now stale
  _svdsComputed = false;
}

#if 0
///////////////////////////////////////////////////////////////////////
// get all the velocity gradients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFdots(const VECTOR& velocity)
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Fs.size() == _tets.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _tets.size(); x++)
  {
    const VECTOR4I& tet = _tets[x];
    VECTOR3 v[4];
    for (int y = 0; y < 4; y++)
    {
      v[y][0] = velocity[3 * tet[y]];
      v[y][1] = velocity[3 * tet[y] + 1];
      v[y][2] = velocity[3 * tet[y] + 2];
    }

    MATRIX3 V;
    V.col(0) = v[1] - v[0];
    V.col(1) = v[2] - v[0];
    V.col(2) = v[3] - v[0];
    _Fdots[x] = V * _DmInvs[x];
  }
}
#endif

///////////////////////////////////////////////////////////////////////
// get the SVD of all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSVDs()
{
  assert(_Us.size() == _triangles.size());
  assert(_Sigmas.size() == _triangles.size());
  assert(_Vs.size() == _triangles.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _triangles.size(); x++)
    svd(_Fs[x], _Us[x], _Sigmas[x], _Vs[x]);

  _svdsComputed = true;
}

///////////////////////////////////////////////////////////////////////
// get deformation gradient
///////////////////////////////////////////////////////////////////////
MATRIX3x2 TRIANGLE_MESH::computeF(const int index) const
{
  const VECTOR3I& t = _triangles[index];
  MATRIX3x2 Ds;
  Ds.col(0) = _vertices[t[1]] - _vertices[t[0]];
  Ds.col(1) = _vertices[t[2]] - _vertices[t[0]];

  return Ds * _DmInvs[index];
}

///////////////////////////////////////////////////////////////////////
// get stretching energy over entire mesh
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeStretchingEnergy(const SHELL::STRETCHING& stretching) const
{
  assert(_triangles.size() == _restTriangleAreas.size());

  VECTOR triangleEnergies(_triangles.size());
  for (int index = 0; index < int(_triangles.size()); index++)
  {
    const MATRIX3x2 F = computeF(index);
    triangleEnergies[index] = _restTriangleAreas[index] * stretching.psi(F);
  }

  return triangleEnergies.sum();
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeStretchingForces(const SHELL::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR9> perElementForces(_triangles.size());
  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const MATRIX3x2& F = _Fs[index];
    const MATRIX3x2 PK1 = stretching.PK1(F);
    const VECTOR9 forceDensity = _pFpxs[index].transpose() * flatten(PK1);
    const VECTOR9 force = -_restTriangleAreas[index] * forceDensity;
    perElementForces[index] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const VECTOR3I& triangle = _triangles[index];
    const VECTOR9& triangleForce = perElementForces[index];
    for (int x = 0; x < 3; x++)
    {
      unsigned int index = 3 * triangle[x];
      forces[index]     += triangleForce[3 * x];
      forces[index + 1] += triangleForce[3 * x + 1];
      forces[index + 2] += triangleForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// get bending energy over entire mesh
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeBendingEnergy(const SHELL::BENDING& bending) const
{
  assert(_triangles.size() == _restTriangleAreas.size());

  VECTOR flapEnergies(_flaps.size());
  for (int index = 0; index < int(_flaps.size()); index++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[index][j]]);

    const REAL theta = _restThetas[index];
    flapEnergies[index] = _restFlapAreas[index] * bending.psi(flap, theta);
  }

  return flapEnergies.sum();
}

///////////////////////////////////////////////////////////////////////
// Use the material gradient to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeBendingForces(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perFlapForces(_flaps.size());
  for (unsigned int index = 0; index < _flaps.size(); index++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[index][j]]);

    const REAL theta = _restThetas[index];
    const VECTOR12 force = -_restFlapAreas[index] * bending.gradient(flap, theta);
    perFlapForces[index] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap      = _flaps[i];
    const VECTOR12& flapForce = perFlapForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * flap[x];

      if (index + 2 >= DOFs)
      {
        cout << " flap index: " << i << endl;
        cout << " vertex index; " << index << endl;
        cout << " DOFs:  " << DOFs << endl;
      }
      assert(index + 2 < DOFs);
      forces[index]     += flapForce[3 * x];
      forces[index + 1] += flapForce[3 * x + 1];
      forces[index + 2] += flapForce[3 * x + 2];
    }
  }
  
  return forces;
}

#if 0
///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the damping force
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeDampingForces(const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const MATRIX3& F = _Fs[tetIndex];
    const MATRIX3& Fdot = _Fdots[tetIndex];
    const MATRIX3 PK1 = damping.PK1(F, Fdot);
    const VECTOR12 forceDensity = _pFpxs[tetIndex].transpose() * flatten(PK1);
    const VECTOR12 force = -_restTetVolumes[tetIndex] * forceDensity;
    perElementForces[tetIndex] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const VECTOR4I& tet = _tets[tetIndex];
    const VECTOR12& tetForce = perElementForces[tetIndex];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the elastic and damping forces
// at the same time
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeInternalForces(const VOLUME::HYPERELASTIC& hyperelastic,
                                       const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perElementForces(_tets.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const MATRIX3& U       = _Us[tetIndex];
    const MATRIX3& V       = _Vs[tetIndex];
    const VECTOR3& Sigma   = _Sigmas[tetIndex];
    const MATRIX3& F = _Fs[tetIndex];
    const MATRIX3& Fdot = _Fdots[tetIndex];

    const MATRIX3 elasticPK1 = hyperelastic.PK1(U, Sigma, V);
    const MATRIX3 dampingPK1 = damping.PK1(F, Fdot);
    const VECTOR12 forceDensity = _pFpxs[tetIndex].transpose() * flatten(elasticPK1 + dampingPK1);
    const VECTOR12 force = -_restTetVolumes[tetIndex] * forceDensity;
    perElementForces[tetIndex] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int tetIndex = 0; tetIndex < _tets.size(); tetIndex++)
  {
    const VECTOR4I& tet = _tets[tetIndex];
    const VECTOR12& tetForce = perElementForces[tetIndex];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the damping gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeDampingHessian(const VOLUME::DAMPING& damping) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& F       = _Fs[i];
    const MATRIX3& Fdot    = _Fdots[i];
    const MATRIX9x12& pFpx = _pFpxs[i];
    const MATRIX9 hessian  = -_restTetVolumes[i] * damping.hessian(F, Fdot);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}
#endif

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeBendingHessian(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perFlapHessians(_flaps.size());
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[i][j]]);

    const REAL theta = _restThetas[i];
    const MATRIX12 hessian  = -_restFlapAreas[i] * bending.hessian(flap, theta);
    perFlapHessians[i] = hessian;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap = _flaps[i];
    const MATRIX12& H = perFlapHessians[i];
    for (int y = 0; y < 4; y++)
    {
      const int yVertex = flap[y];
      for (int x = 0; x < 4; x++)
      {
        const int xVertex = flap[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeBendingClampedHessian(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perFlapHessians(_flaps.size());
//#pragma omp parallel
//#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[i][j]]);

    const REAL theta = _restThetas[i];
    const MATRIX12 hessian  = -_restFlapAreas[i] * bending.clampedHessian(flap, theta);
    perFlapHessians[i] = hessian;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap = _flaps[i];
    const MATRIX12& H = perFlapHessians[i];
    for (int y = 0; y < 4; y++)
    {
      const int yVertex = flap[y];
      for (int x = 0; x < 4; x++)
      {
        const int xVertex = flap[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeStretchingHessian(const SHELL::STRETCHING& stretching) const
{
  vector<MATRIX9> perElementHessians(_triangles.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const MATRIX3x2& F       = _Fs[i];
    const MATRIX6x9& pFpx = _pFpxs[i];
    const MATRIX6 hessian  = -_restTriangleAreas[i] * stretching.hessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I& triangle = _triangles[i];
    const MATRIX9& H = perElementHessians[i];
    for (int y = 0; y < 3; y++)
    {
      const int yVertex = triangle[y];
      for (int x = 0; x < 3; x++)
      {
        const int xVertex = triangle[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeStretchingClampedHessian(const SHELL::STRETCHING& stretching) const
{
  vector<MATRIX9> perElementHessians(_triangles.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const MATRIX3x2& F       = _Fs[i];
    const MATRIX6x9& pFpx = _pFpxs[i];
    const MATRIX6 hessian  = -_restTriangleAreas[i] * stretching.clampedHessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I& triangle = _triangles[i];
    const MATRIX9& H = perElementHessians[i];
    for (int y = 0; y < 3; y++)
    {
      const int yVertex = triangle[y];
      for (int x = 0; x < 3; x++)
      {
        const int xVertex = triangle[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setPositions(const VECTOR& positions)
{
  assert(positions.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = positions[3 * x];
    _vertices[x][1] = positions[3 * x + 1];
    _vertices[x][2] = positions[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// get the current displacement in vector form
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::getDisplacement() const
{
  VECTOR delta(_vertices.size() * 3);
  delta.setZero();

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    const VECTOR3 diff = _vertices[x] - _restVertices[x];
    const int x3 = 3 * x;
    delta[x3]     = diff[0];
    delta[x3 + 1] = diff[1];
    delta[x3 + 2] = diff[2];
  }

  return delta;
}

///////////////////////////////////////////////////////////////////////
// set the vertex displacements to these values exactly
// ALSO UPDDATES _verticesOld
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3); 
  std::swap(_vertices,_verticesOld); //Update the old vertices 
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = _restVertices[x][0] + delta[3 * x];
    _vertices[x][1] = _restVertices[x][1] + delta[3 * x + 1];
    _vertices[x][2] = _restVertices[x][2] + delta[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// add the following deltas to the position
// DOES NOT UPDATE _verticesOld
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] += delta[3 * x];
    _vertices[x][1] += delta[3 * x + 1];
    _vertices[x][2] += delta[3 * x + 2];
  }
}

#if 0
///////////////////////////////////////////////////////////////////////
// get the bounding box for the current mesh
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::getBoundingBox(VECTOR3& mins, VECTOR3& maxs) const
{
  assert(_vertices.size() > 0);
  mins = _vertices[0];
  maxs = _vertices[0];

  for (unsigned int x = 1; x < _vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (mins[y] < _vertices[x][y]) ? mins[y] : _vertices[x][y];
      maxs[y] = (maxs[y] > _vertices[x][y]) ? maxs[y] : _vertices[x][y];
    }
}

///////////////////////////////////////////////////////////////////////
// write out the surface to OBJ triangle mesh
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeSurfaceToObj(const string& filename, const TRIANGLE_MESH& tetMesh)
{
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << " Failed to open file " << filename.c_str() << endl;
    return false;
  }
  cout << " Writing out file " << filename.c_str() << " ... " << flush;

  const vector<VECTOR3>& vertices = tetMesh.vertices();
  const vector<VECTOR3I>& surfaceTriangles = tetMesh.surfaceTriangles();

  // do the ugly thing and just write out all the vertices, even
  // the internal ones
  for (unsigned int x = 0; x < vertices.size(); x++)
    fprintf(file, "v %f %f %f\n", vertices[x][0], vertices[x][1], vertices[x][2]);

  // write out the indices for the surface triangles, but remember
  // that Objs are 1-indexed
  for (unsigned int x = 0; x < surfaceTriangles.size(); x++)
    fprintf(file, "f %i %i %i\n", surfaceTriangles[x][0] + 1, 
                                  surfaceTriangles[x][1] + 1,
                                  surfaceTriangles[x][2] + 1);

  fclose(file);
  cout << " done. " << endl;
  return true;
}
#endif

///////////////////////////////////////////////////////////////////////
// read in OBJ-style tet mesh file using Tiny OBJ Loader
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readObjFile(const string& filename, 
                                vector<VECTOR3>& vertices,
                                vector<VECTOR3I>& triangles)
{
  // load up the OBJ
  cout << " Reading in *.obj file " << filename.c_str() << endl;

  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  std::string warn;
  std::string err;
  const bool triangulate = true;
  bool success = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str(),
                                  NULL, triangulate);

  if (!warn.empty()) {
    cout << "WARNING: " << warn << endl;
  }
  if (!err.empty()) {
    cout << "ERROR: " << err << endl;
  }
  if (!success) {
    cout << " Failed to read " << filename.c_str() << endl;
    return false;
  }

  // erase whatever was in the vectors before
  vertices.clear();
  triangles.clear();

  // store up the vertices
  const vector<float>& vs = attrib.vertices;
  for (unsigned int x = 0; x < (vs.size() / 3); x++)
  {
    const int x3 = 3 * x;
    const VECTOR3 v(vs[x3 + 0], vs[x3 + 1], vs[x3 + 2]);
    vertices.push_back(v);
  }
  cout << " Read in " << vs.size() / 3 << " vertices " << endl;

  // there's at least one shape, right?
  assert(shapes.size() > 0);

  // just pick off the first shape
  size_t index_offset = 0;
  for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++)
  {
    VECTOR3I triangle;

    // it's a triangle, right?
    const size_t fnum = shapes[0].mesh.num_face_vertices[f];
    assert(fnum == 3);
    for (size_t v = 0; v < 3; v++) 
    {
      tinyobj::index_t idx = shapes[0].mesh.indices[index_offset + v];
      int index = idx.vertex_index;

      triangle[v] = index;
    }
    index_offset += 3;

    triangles.push_back(triangle);
  }
  cout << " Read in " << triangles.size() << " faces " << endl;

  return true;
}


///////////////////////////////////////////////////////////////////////
// write out OBJ-style tet mesh file
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeObjFile(const string& filename, const TRIANGLE_MESH& triMesh, 
                             const bool restVertices)
{
  FILE* file = fopen(filename.c_str(), "w");
  cout << " Writing out mesh file: " << filename.c_str() << endl;
 
  if (file == NULL)
  {
    cout << " Failed to open file!" << endl;
    return false;
  }

  const vector<VECTOR3>& vertices = (restVertices) ? triMesh.restVertices() : triMesh.vertices(); 
  const vector<VECTOR3I>& tris = triMesh.triangles(); 

  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    const VECTOR3& v = vertices[x];
    fprintf(file, "v %.17g %.17g %.17g\n", v[0], v[1], v[2]);
  }
  for (unsigned int x = 0; x < tris.size(); x++)
  {
    const VECTOR3I& tri = tris[x];
    fprintf(file, "t %i %i %i %i\n", tri[0], tri[1], tri[2]);
  }

  fclose(file);
  return true;
}

///////////////////////////////////////////////////////////////////////
// normalize vertices so that they're in a unit box, 
// centered at (0.5, 0.5, 0.5)
///////////////////////////////////////////////////////////////////////
vector<VECTOR3> TRIANGLE_MESH::normalizeVertices(const vector<VECTOR3>& vertices)
{
  assert(vertices.size() > 0);
  VECTOR3 mins = vertices[0];
  VECTOR3 maxs = vertices[0];
  for (unsigned int x = 1; x < vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (mins[y] < vertices[x][y]) ? mins[y] : vertices[x][y];
      maxs[y] = (maxs[y] > vertices[x][y]) ? maxs[y] : vertices[x][y];
    }

  const VECTOR3 lengths = maxs - mins;
  const REAL maxLengthInv = 1.0 / lengths.maxCoeff();

  vector<VECTOR3> normalized = vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    normalized[x] -= mins;
    normalized[x] *= maxLengthInv;
    
    normalized[x] += VECTOR3(0.5, 0.5, 0.5);
  }

  return normalized;
}

///////////////////////////////////////////////////////////////////////
// see if the projection of v onto the plane of v0,v1,v2 is inside 
// the triangle formed by v0,v1,v2
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                           const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute distance between a point and triangle
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                     const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
  {
    const VECTOR3 nHat = n / n.norm();
    const REAL normalDistance = (nHat.dot(v - v0));
    return fabs(normalDistance);
  }

  // project onto each edge, find the distance to each edge
  const VECTOR3 e3 = v2 - v1;
  const VECTOR3 ev = v - v0;
  const VECTOR3 ev3 = v - v1;
  const VECTOR3 e1Hat = e1 / e1.norm();
  const VECTOR3 e2Hat = e2 / e2.norm();
  const VECTOR3 e3Hat = e3 / e3.norm();
  VECTOR3 edgeDistances(FLT_MAX, FLT_MAX, FLT_MAX);

  // see if it projects onto the interval of the edge
  // if it doesn't, then the vertex distance will be smaller,
  // so we can skip computing anything
  const REAL e1dot = e1Hat.dot(ev);
  if (e1dot > 0.0 && e1dot < e1.norm())
  {
    const VECTOR3 projected = v0 + e1Hat * e1dot;
    edgeDistances[0] = (v - projected).norm();
  }
  const REAL e2dot = e2Hat.dot(ev);
  if (e2dot > 0.0 && e2dot < e2.norm())
  {
    const VECTOR3 projected = v0 + e2Hat * e2dot;
    edgeDistances[1] = (v - projected).norm();
  }
  const REAL e3dot = e3Hat.dot(ev3);
  if (e3dot > 0.0 && e3dot < e3.norm())
  {
    const VECTOR3 projected = v1 + e3Hat * e3dot;
    edgeDistances[2] = (v - projected).norm();
  }

  // get the distance to each vertex
  const VECTOR3 vertexDistances((v - v0).norm(), 
                                (v - v1).norm(), 
                                (v - v2).norm());

  // get the smallest of both the edge and vertex distances
  const REAL vertexMin = vertexDistances.minCoeff();
  const REAL edgeMin = edgeDistances.minCoeff();

  // return the smallest of those
  return (vertexMin < edgeMin) ? vertexMin : edgeMin;
}

///////////////////////////////////////////////////////////////////////
// see if the vertex is inside the collision cell described in 
// Chapter 11: Collision Processing, [Kim and Eberle 2020]
//
// this seems to work okay, but need stronger debugging checks
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::insideCollisionCell(const int triangleID, const VECTOR3& vertex)
{
  const VECTOR3I& t = _triangles[triangleID];
  vector<VECTOR3> v;
  v.push_back(_vertices[t[0]]);
  v.push_back(_vertices[t[1]]);
  v.push_back(_vertices[t[2]]);
  VECTOR3 n = planeNormal(v);

  // get the normals of the three adjacent faces
  const VECTOR3I& neighbors = _triangleNeighbors[triangleID];
  for (int x = 0; x < 3; x++)
  {
    VECTOR3 nebHat;
    // Just skip if it's boundary edge; distance checks in higher function
    // will make sure nothing crazy happens
    if (neighbors[x] < 0){
      // int j = -neighbors[x] -1;
      // VECTOR3 p0 = (v[(j+1)%3] + v[j])/2;
      // VECTOR3 p1 = v[(j+2)%3];
      // VECTOR3 edge = (v[j+1%3] - v[j]).normalized();
      // nebHat = (p1-p0 -(p1-p0).dot(edge)*edge).normalized();
      continue;
    }
    else{
      const VECTOR3I& tNeighbor = _triangles[neighbors[x]];
      vector<VECTOR3> vNeighbor;
      vNeighbor.push_back(_vertices[tNeighbor[0]]);
      vNeighbor.push_back(_vertices[tNeighbor[1]]);
      vNeighbor.push_back(_vertices[tNeighbor[2]]);
      const VECTOR3 nNeighbor = planeNormal(vNeighbor);
      // do the inside check
      const VECTOR3 ne = (nNeighbor + n).normalized();
      const VECTOR3 eij = v[(x+1) % 3] - v[x];
      const VECTOR3 neb = ne.cross(eij);
      nebHat = neb.normalized();
    }

    const REAL deplane = nebHat.dot(vertex - v[x]);
    if (deplane < 0.0)
      return false;
  }

  // TODO: some face plane compares, and returning if the vertex is above
  // or below the face
  return true;
}

// ///////////////////////////////////////////////////////////////////////
// // compute distance to collision cell wall, where positive means inside
// // and negative means outside
// // does this get used anywhere?
// ///////////////////////////////////////////////////////////////////////
// REAL TRIANGLE_MESH::distanceToCollisionCellWall(const int triangleID, const VECTOR3& vertex)
// {
//   const VECTOR3I& t = _triangles[triangleID];
//   vector<VECTOR3> v;
//   v.push_back(_vertices[t[0]]);
//   v.push_back(_vertices[t[1]]);
//   v.push_back(_vertices[t[2]]);
//   VECTOR3 n = planeNormal(v);

//   REAL smallestDistance = FLT_MAX;
//   // get the normals of the three adjacent faces
//   vector<VECTOR3> nNeighbors;
//   const VECTOR3I& neighbors = _triangleNeighbors[triangleID];
//   for (int x = 0; x < 3; x++)
//   {
//     VECTOR3 nebHat;
//     // inside check is pretty different if considering orphaned or unorphaned edge
//     if (neighbors[x] < 0){
//       int j = -neighbors[x] -1;
//       VECTOR3 p0 = (v[(j+1)%3] + v[j])/2;
//       VECTOR3 p1 = v[(j+2)%3];
//       nebHat = (p1-p0).normalized();
//     }
//     else{
//       const VECTOR3I& tNeighbor = _triangles[neighbors[x]];
//       vector<VECTOR3> vNeighbor;
//       vNeighbor.push_back(_vertices[tNeighbor[0]]);
//       vNeighbor.push_back(_vertices[tNeighbor[1]]);
//       vNeighbor.push_back(_vertices[tNeighbor[2]]);
//       const VECTOR3 nNeighbor = planeNormal(vNeighbor);
//       // do the inside check
//       const VECTOR3 ne = (nNeighbor + n).normalized();
//       const VECTOR3 eij = v[(x+1) % 3] - v[x];
//       const VECTOR3 neb = ne.cross(eij);
//       nebHat = neb.normalized();
//     }

//     const REAL deplane = nebHat.dot(vertex - v[x]);
//     if (fabs(deplane) < smallestDistance)
//     {
//       smallestDistance = fabs(deplane);
//     }
//   }

//   return smallestDistance;
// }


///////////////////////////////////////////////////////////////////////
// find all the vertex-face collision pairs, using the 
// InFaceRegion test from "Collision Processing" chapter of
// "Dynamic Deformables"
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexFaceCollisions()
{
  TIMER functionTimer(__FUNCTION__);

  // if a vertex is part of an inverted triangle, don't have it participate 
  // in a self-collision. That triangle needs to get its house in order 
  // before it starts bossing around a surface face. Not checking for 
  // this causes faces to get horribly tangled in inverted configurations.
  computeInvertedVertices();

  _vertexFaceCollisions.clear();
  const REAL collisionEps = _collisionEps;

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    // if the vertex is involved in an inverted triangle, give up
    if (_invertedVertices[x]) 
      continue;

    const VECTOR3& surfaceVertex = _vertices[x];

    // find the close triangles
    for (unsigned int y = 0; y < _triangles.size(); y++)
    {
      // if the surface triangle is so small the normal could be degenerate, skip it
      if (surfaceTriangleIsDegenerate(y))
        continue;

      const VECTOR3I& t = _triangles[y];

      // if it's an inverted face, move on
      if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
        continue;

      // if this triangle is in the one-ring of the current vertex, skip it
      if (t[0] == x || t[1] == x || t[2] == x) continue;
      
      const REAL distance = pointTriangleDistance(_vertices[t[0]], _vertices[t[1]],
                                                  _vertices[t[2]], surfaceVertex);

      if (distance < collisionEps)
      {
        // if the point, projected onto the face's plane, is inside the face,
        // then record the collision now
        if (pointProjectsInsideTriangle(_vertices[t[0]], _vertices[t[1]], 
                                        _vertices[t[2]], surfaceVertex))
        {
          pair<int,int> collision(x, y);
          _vertexFaceCollisions.push_back(collision);
          continue;
        }
        if (insideCollisionCell(y, surfaceVertex))
        {
          pair<int,int> collision(x, y);
          _vertexFaceCollisions.push_back(collision);
        }
      }
    }
  }

#if VERY_VERBOSE
  if (_vertexFaceCollisions.size() > 0){
    cout << " Found " << _vertexFaceCollisions.size() << " vertex-face collisions " << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const pair<int,int> collision = _vertexFaceCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// find vertex-face CCD collisions 
// through comparing _verticesOld with _vertices
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexFaceCCD()
{
#ifndef DISABLE_CCD
  TIMER functionTimer(__FUNCTION__);

  // if a vertex is part of an inverted triangle, don't have it participate 
  // in a self-collision. That triangle needs to get its house in order 
  // before it starts bossing around a surface face. Not checking for 
  // this causes faces to get horribly tangled in inverted configurations.
  computeInvertedVertices();

  _vertexFaceCCDs.clear(); // data structure for storing CCD vf collisions
  const REAL collisionEps = _collisionEps;

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    // if the vertex is involved in an inverted triangle, give up
    if (_invertedVertices[x]) 
      continue;

    // find the close triangles
    for (unsigned int y = 0; y < _triangles.size(); y++)
    {
      // if the surface triangle is so small the normal could be degenerate, skip it
      if (surfaceTriangleIsDegenerate(y))
        continue;

      const VECTOR3I& t = _triangles[y];

      // if it's an inverted face, move on
      if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
        continue;

      // if this triangle is in the one-ring of the current vertex, skip it
      if (t[0] == x || t[1] == x || t[2] == x) continue;
      
      // Do the CCD test: If CCD succeeds, save the vertex index and the triangle index
      
      // The colliding vertex
      REAL* x01 = _vertices[x].data();
      REAL* x00 = _verticesOld[x].data();
      // The triangle
      REAL* x10 = _verticesOld[t[0]].data();
      REAL* x11 = _vertices[t[0]].data();
      REAL* x20 = _verticesOld[t[1]].data();
      REAL* x21 = _vertices[t[1]].data();
      REAL* x30 = _verticesOld[t[2]].data();
      REAL* x31 = _vertices[t[2]].data();

      REAL tHit;

      if(_ccd->Vertex_Triangle_CCD(x00,x01,x10,x11,x20,x21,x30,x31,tHit)){
        pair<REAL,pair<int,int>> collisionPair;
        collisionPair.second.first = x;
        collisionPair.second.second = y;
        collisionPair.first = tHit;
        _vertexFaceCCDs.push_back(collisionPair);
      }

    }
  }

  // sort the vf CCDs by time
  std::sort(_vertexFaceCCDs.begin(), _vertexFaceCCDs.end());

#if VERY_VERBOSE
  if (_vertexFaceCCDs.size() > 0){
    cout << " Found " << _vertexFaceCCDs.size() << " vertex-face CCDs " << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _vertexFaceCCDs.size(); x++)
  {
    const pair<REAL, pair<int,int>> collision = _vertexFaceCCDs[x];
    cout << "(" << collision.first << ", " << collision.second.first << ", " << collision.second.second << ")" << endl;
  }
#endif
  #endif
  return;
}

///////////////////////////////////////////////////////////////////////
// find all the edge-edge self collision pairs, using the 
// brute-force tests
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeEdgeCollisions()
{
  TIMER functionTimer(__FUNCTION__);
  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  //_edgeEdgeCollisionEps.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();

  // build a mapping from edge index pairs to _edges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    edgeHash[edge] = x;
  }

  // get the nearest edge to each edge, not including itself
  // and ones where it shares a vertex
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    int closestEdge = -1;
    REAL closestDistance = FLT_MAX;
    VECTOR2 aClosest(-1,-1);
    VECTOR2 bClosest(-1,-1);
    const VECTOR2I outerEdge = _edges[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];

    //const unsigned int outerFlat = outerEdge[0] + outerEdge[1] * _edges.size();
    
    // find the closest other edge
    for (unsigned int y = x + 1; y < _edges.size(); y++)
    {
      const VECTOR2I innerEdge = _edges[y];
      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;

      const VECTOR3& v2 = _vertices[innerEdge[0]];
      const VECTOR3& v3 = _vertices[innerEdge[1]];

      /*
      // call the geometric test
      VECTOR3 innerPoint, outerPoint, midpoint, normal;
      bool intersect;
      IntersectLineSegments(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2],
                            v2[0], v2[1], v2[2], v3[0], v3[1], v3[2],
                            false, 1e-8,
                            outerPoint[0], outerPoint[1], outerPoint[2],
                            innerPoint[0], innerPoint[1], innerPoint[2],
                            midpoint[0], midpoint[1], midpoint[2],
                            normal[0], normal[1], normal[2], intersect);
                            */
      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);  

      const REAL distance = (innerPoint - outerPoint).norm();

      //if (distance < separationDistance)
      //  separationDistance = distance;

      if (distance > closestDistance) continue;

      // get the line interpolation coordinates
      VECTOR2 a,b;
      const VECTOR3 e0 = v1 - v0;
      const VECTOR3 e1 = v3 - v2;

      // this is a little dicey in general, but if the intersection test isn't
      // total garbage, it should still be robust
      a[1] = (outerPoint - v0).norm() / e0.norm();
      a[0] = 1.0 - a[1];
      b[1] = (innerPoint - v2).norm() / e1.norm();
      b[0] = 1.0 - b[1];

      // if it's really close to an end vertex, skip it
      const REAL skipEps = 1e-4;
      if ((a[0] < skipEps) || (a[0] > 1.0 - skipEps)) continue;
      if ((a[1] < skipEps) || (a[1] > 1.0 - skipEps)) continue;
      if ((b[0] < skipEps) || (b[0] > 1.0 - skipEps)) continue;
      if ((b[1] < skipEps) || (b[1] > 1.0 - skipEps)) continue;

      // it's mid-segment, and closest, so remember it
      closestDistance = distance;
      closestEdge = y;

      aClosest = a;
      bClosest = b;
    }

    // if nothing was close, move on
    if (closestEdge == -1) continue;

    /*
    // retrieve the eps of the closest edge
    const VECTOR2I innerEdge = _edges[closestEdge];
    const unsigned int innerFlat = innerEdge[0] + innerEdge[1] * _edges.size();
    const pair<unsigned int, unsigned int> edgeEdge(innerFlat, outerFlat);

    // it exists, right?
    assert(_edgeEdgeRestDistance.find(edgeEdge) != _edgeEdgeRestDistance.end());
    //const REAL collisionEps = _edgeEdgeRestDistance[edgeEdge];
    const REAL eeDistance = _edgeEdgeRestDistance[edgeEdge];
    const REAL collisionEps = (eeDistance < _collisionEps) ? eeDistance : _collisionEps;
    */

    // are they within each other's one rings?
    const VECTOR2I innerEdge = _edges[closestEdge];
    bool insideOneRing = false;

    for (int j = 0; j < 2; j++)
    {
      pair<int, int> lookup;
      lookup.first = outerEdge[j];
      for (int i = 0; i < 2; i++)
      {
        lookup.second = innerEdge[i];
        if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end())
          insideOneRing = true;
      }
    }
    if (insideOneRing) continue;

    /*
    const VECTOR2I innerEdge = _edges[closestEdge];
    pair<int, int> lookup;
    lookup.first = outerEdge[0];
    lookup.second = innerEdge[0];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[0];
    lookup.second = innerEdge[1];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[1];
    lookup.second = innerEdge[0];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[1];
    lookup.second = innerEdge[1];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    */

    // if it's within the positive threshold, it's in collision
    if (closestDistance < _collisionEps)
    {
      pair<int,int> collision(x, closestEdge);
      _edgeEdgeCollisions.push_back(collision);
      
      // this was actually set, right?
      assert(aClosest[0] > 0.0 && aClosest[1] > 0.0);
      assert(bClosest[0] > 0.0 && bClosest[1] > 0.0);

      pair<VECTOR2,VECTOR2> coordinate(aClosest, bClosest);
      _edgeEdgeCoordinates.push_back(coordinate);

      // get the areas too
      const VECTOR2I innerEdge = _edges[closestEdge];
      const pair<int,int> outerPair(outerEdge[0], outerEdge[1]);
      const pair<int,int> innerPair(innerEdge[0], innerEdge[1]);
      const REAL xArea = _restEdgeAreas[edgeHash[outerPair]];
      const REAL closestArea = _restEdgeAreas[edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);

      // find out if they are penetrating
      vector<VECTOR3> edge(2);
      edge[0] = v0;
      edge[1] = v1;

      // get the adjacent triangles of the *other* edge
      VECTOR2I adjacentTriangles = _edgeTriangleNeighbors[edgeHash[innerPair]];

      // build triangle 0
      const VECTOR3I surfaceTriangle0 = _triangles[adjacentTriangles[0]];
      vector<VECTOR3> triangle0;
      triangle0.push_back(_vertices[surfaceTriangle0[0]]);
      triangle0.push_back(_vertices[surfaceTriangle0[1]]);
      triangle0.push_back(_vertices[surfaceTriangle0[2]]);

      // build triangle 1
      vector<VECTOR3> triangle1;
      // if there's another triangle on the other side (this is in case we're looking at cloth)
      // then store that one too
      if (adjacentTriangles[1] != -1)
      {
        const VECTOR3I surfaceTriangle1 = _triangles[adjacentTriangles[1]];
        triangle1.push_back(_vertices[surfaceTriangle1[0]]);
        triangle1.push_back(_vertices[surfaceTriangle1[1]]);
        triangle1.push_back(_vertices[surfaceTriangle1[2]]);
      }

      // see if the edges are already penetrating the opposing faces
      bool penetrating = false;
      if (triangle0.size() > 0) penetrating = faceEdgeIntersection(triangle0, edge);
      if (triangle1.size() > 0) penetrating = penetrating || faceEdgeIntersection(triangle1, edge);

      _edgeEdgeIntersections.push_back(penetrating);
      //_edgeEdgeCollisionEps.push_back(_collisionEps);

      // TODO: for completeness, should probably test the other edges against the other
      // pair, just in case we're looking at a degenerate case. In general, seems redundant.
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());

#if VERY_VERBOSE
  if (_edgeEdgeCollisions.size() > 0){
    cout << " Collision area array size: " << _edgeEdgeCollisionAreas.size() << endl;
    cout << " Intersections array size: " << _edgeEdgeIntersections.size() << endl;
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _edgeEdgeCollisions.size(); x++)
  {
    const pair<int,int> collision = _edgeEdgeCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// find all the edge-edge CCD collisions
// brute-force tests 
// by comparison of _verticesOld and _vertices
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeEdgeCCD()
{
  #ifndef DISABLE_CCD
  TIMER functionTimer(__FUNCTION__);
  _edgeEdgeCCDs.clear();

  // build a mapping from edge index pairs to _edges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    edgeHash[edge] = x;
  }

  // dumb loop through all the edges and manually CCD-ing
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    VECTOR2I outerEdge = _edges[x];
    for (unsigned int y = x + 1; y < _edges.size(); y++)
    {
      VECTOR2I innerEdge = _edges[y];
      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;
      //Do edge-edge CCD and save to _edgeEdgeCCDs
      //Outer edge
      REAL* x01 = _vertices[outerEdge[0]].data();
      REAL* x00 = _verticesOld[outerEdge[0]].data();
      REAL* x10 = _verticesOld[outerEdge[1]].data();
      REAL* x11 = _vertices[outerEdge[1]].data();
      //Inner edge
      REAL* x20 = _verticesOld[innerEdge[0]].data();
      REAL* x21 = _vertices[innerEdge[0]].data();
      REAL* x30 = _verticesOld[innerEdge[1]].data();
      REAL* x31 = _vertices[innerEdge[1]].data();

      REAL tHit;

      if(_ccd->Edge_Edge_CCD(x00,x01,x10,x11,x20,x21,x30,x31,tHit)){
        pair< REAL,pair<int,int>> collision;
        collision.second.first = x;
        collision.second.second = y;
        collision.first = tHit;
        _edgeEdgeCCDs.push_back(collision);
      }

    }

  }

  // sort the ee CCDs by time
  std::sort(_edgeEdgeCCDs.begin(),_edgeEdgeCCDs.end());

#if VERY_VERBOSE
  if (_edgeEdgeCCDs.size() > 0){
    cout << " Found " << _edgeEdgeCCDs.size() << " edge-edge CCDs " << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _edgeEdgeCCDs.size(); x++)
  {
    const pair<REAL, pair<int,int>> collision = _edgeEdgeCCDs[x];
    cout << "(" << collision.first << ", " << collision.second.first << ", " << collision.second.second <<  ")" << endl;
  }
#endif
  #endif
  return;
}

///////////////////////////////////////////////////////////////////////
// gather and sort all ee and vf CCDs into a single vector
// first is 0 for vf, 1 for ee; second is VECTOR4I of vertices involved
// vf ordering: ensured such that first vert is the point, last three
// oriented triangle st normal points outwards (according to verticesold)
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::gatherCCDs()
{
  TIMER functionTimer(__FUNCTION__);
  // Merge sorting very cumbersome data structures with extra metadata
  _allCCDs.clear();
  int mvf = 0;
  int mee = 0;
  bool vfDepleted = mvf >= (_vertexFaceCCDs.size());
  bool eeDepleted = mee >= (_edgeEdgeCCDs.size());

  while(!vfDepleted || !eeDepleted){
    // Initializing vfTime and eeTime >1 so 
    // other time will always be lower
    // (catching if one vector depletes faster than the other)
    pair<REAL,VECTOR5I> ccdee; 
    pair<REAL,VECTOR5I> ccdvf;
    REAL vfTime = 2.0;
    REAL eeTime = 2.0;
    if(!eeDepleted){ //construct the ee pair (if exists)
      eeTime = _edgeEdgeCCDs[mee].first;
      ccdee.first = eeTime;
      VECTOR2I e0 = _edges[_edgeEdgeCCDs[mee].second.first];
      VECTOR2I e1 = _edges[_edgeEdgeCCDs[mee].second.second];
      ccdee.second[0] = e0[0];
      ccdee.second[1] = e0[1];
      ccdee.second[2] = e1[0];
      ccdee.second[3] = e1[1];
      ccdee.second[4] = 1;
    }
    if(!vfDepleted){ //construct the vf pair (if exists)
      vfTime = _vertexFaceCCDs[mvf].first;
      ccdvf.first = vfTime;
      const int vertexID = _vertexFaceCCDs[mvf].second.first;
      const int faceID = _vertexFaceCCDs[mvf].second.second;
      const VECTOR3I& face = _triangles[faceID];
    
      vector<VECTOR3> vs(4);
      for(int i = 0; i < 3; i++){
        vs[i] = _verticesOld[face[i]];
      }

      // build a tet with the correct vertex ordering
      VECTOR4I tet;
      tet[0] = vertexID;

      // cloth collisions: if face normal is aligned with collision
      // edge, invert face order. Otherwise feed face order normally

      VECTOR3 faceNorm = planeNormal(vs);
      VECTOR3 wTests = _verticesOld[vertexID] - _verticesOld[face[1]];
      if(wTests.dot(faceNorm) > 0)
      {
        tet[1] = face[2];
        tet[2] = face[1];
        tet[3] = face[0];
      }
      else
      {
        tet[1] = face[0];
        tet[2] = face[1];
        tet[3] = face[2];
      }
      ccdvf.second[0] = tet[0];
      ccdvf.second[1] = tet[1];
      ccdvf.second[2] = tet[2];
      ccdvf.second[3] = tet[3];
      ccdvf.second[4] = 0;
    }
    // do the comparison
    if(vfTime <= eeTime){
      _allCCDs.push_back(ccdvf);
      mvf++;
    }
    else{
      _allCCDs.push_back(ccdee);
      mee++;
    }
    vfDepleted = mvf >= (_vertexFaceCCDs.size());
    eeDepleted = mee >= (_edgeEdgeCCDs.size());
  }

  #if VERY_VERBOSE
  if(_allCCDs.size() > 0){
    cout << " Found " << _allCCDs.size() << " total CCDs" << endl;
    cout << "_allCCDs: " << endl;
    for(unsigned int i = 0; i < _allCCDs.size(); i++){
      const pair<REAL,VECTOR5I> ccd = _allCCDs[i];
      cout << "(" << ccd.first << ", " << ccd.second.transpose() << ")" << endl;
    }
  }
  #endif
}

///////////////////////////////////////////////////////////////////////
// gather _allCCDs into regions by checking contiguous vertices
// updates _regionsCCD
// Basically connected component finding but we'll just do something slow
// and easy to implement
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::gatherCCDRegions()
{
  TIMER functionTimer(__FUNCTION__);
  _regionsCCD.clear();
  _vertexRegionsCCD.clear();
  // Loop through _allCCDs once just to get clean _vertexRegionsCCD
  // (Then a second time to assign chronologically to the regions)

  // FIRST LOOP
  for(unsigned int i = 0; i < _allCCDs.size(); i++){
    // The ccd in question
    pair<REAL,VECTOR5I> ccd = _allCCDs[i]; 
    // Accumulator for the new glob this ccd might make
    set<int> accumulator;
    // the set of vertices that the collision represents
    set<int> collVertexSet; 
    collVertexSet.insert(ccd.second[0]);
    collVertexSet.insert(ccd.second[1]);
    collVertexSet.insert(ccd.second[2]);
    collVertexSet.insert(ccd.second[3]);
    vector<bool> staysAlong;
    for(unsigned int j = 0; j < _vertexRegionsCCD.size(); j++){
      // The set of vertices already in the region.
      set<int> regionSet = _vertexRegionsCCD[j];

      // Does the region share any vertices with the ccd vertices?
      set<int> sharedVerts;
      std::set_intersection(regionSet.begin(),regionSet.end(),
                            collVertexSet.begin(),collVertexSet.end(),
                            std::inserter(sharedVerts,sharedVerts.begin()));
      // If so, add its info to the region vertices (but not yet with the region data)
      staysAlong.push_back(true);
      if(sharedVerts.size() !=0){
        accumulator.insert(_vertexRegionsCCD[j].begin(), _vertexRegionsCCD[j].end());
        staysAlong.back() = false;
      }
    }
    if(accumulator.size() == 0){
      _vertexRegionsCCD.push_back(collVertexSet);
    }
    else{
      vector<set<int>> regionReplace;
      accumulator.insert(collVertexSet.begin(),collVertexSet.end());
      regionReplace.push_back(accumulator);
      for(unsigned int k = 0; k < _vertexRegionsCCD.size(); k++){
       if(staysAlong[k]) regionReplace.push_back(_vertexRegionsCCD[k]); 
      }
      _vertexRegionsCCD = regionReplace;
    }
  }

  _regionsCCD.resize(_vertexRegionsCCD.size());

  // SECOND LOOP: All the regions are created and work contiguously; now just add the corresponding CCDs to their regions
  for(unsigned int i = 0; i < _allCCDs.size(); i++){
    // The ccd in question
    pair<REAL,VECTOR5I> ccd = _allCCDs[i]; 
    // the set of vertices that the collision represents
    set<int> collVertexSet; 
    collVertexSet.insert(ccd.second[0]);
    collVertexSet.insert(ccd.second[1]);
    collVertexSet.insert(ccd.second[2]);
    collVertexSet.insert(ccd.second[3]);
    for(unsigned int j = 0; j < _vertexRegionsCCD.size(); j++){
      set<int> regionSet = _vertexRegionsCCD[j];
      // Does the region share any vertices with the ccd vertices?
      set<int> sharedVerts;
      std::set_intersection(regionSet.begin(),regionSet.end(),
                            collVertexSet.begin(),collVertexSet.end(),
                            std::inserter(sharedVerts,sharedVerts.begin()));
      if(sharedVerts.size() !=0){
        _regionsCCD[j].push_back(i);
      }
    }
  }

  #if VERY_VERBOSE
  if(_regionsCCD.size() > 0){
    cout << "Found " << _regionsCCD.size() << "CCD regions" << endl;
    cout << "_regionsCCD.size(): " << _regionsCCD.size() << " | _vertexRegionsCCD.size(): " << _vertexRegionsCCD.size() << endl;
    cout << "_vertexRegionsCCD: " << endl;
    for(const set<int>& vSet : _vertexRegionsCCD){
      cout << "{";
      for(const int& v : vSet){
        cout << v << ", ";
      }
      cout << "}" << endl;
    }
  }
  #endif

}

///////////////////////////////////////////////////////////////////////
// loops through the regions and resolves the first chronologically occurring
// collision
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::solveCCDRegions(){
  for (unsigned int i = 0; i < _regionsCCD.size(); i++){
    //The collision in question 
    const pair<REAL, VECTOR5I> coll = _allCCDs[_regionsCCD[i][0]];
    
    vector<VECTOR3> vs;
    vector<VECTOR3> ve;
    for(int i = 0 ; i < 4; i++){
      ve.push_back(_vertices[coll.second[i]]);
      vs.push_back(_verticesOld[coll.second[i]]);
    }

    REAL tHit = coll.first;
    const int collType = coll.second[4];
    if(collType == 1){ //edge edge UNFINISHED

    }
    else if(collType == 0){ //vertex face UNFINISHED

    }
  }
}

///////////////////////////////////////////////////////////////////////
// get the triangle area
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::triangleArea(const vector<VECTOR3>& triangle)
{
  const VECTOR3 edge1 = triangle[1] - triangle[0];
  const VECTOR3 edge2 = triangle[2] - triangle[0];
  return 0.5 * edge1.cross(edge2).norm();
}


///////////////////////////////////////////////////////////////////////
// get the normal to a plane, specified by three points
// returns normal assuming positive oriented vertex ordering
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
}

///////////////////////////////////////////////////////////////////////
// project point onto plane, specific by three points
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point)
{
  const VECTOR3 normal = planeNormal(plane);
  return point - (normal.dot(point - plane[0])) * normal;
}


///////////////////////////////////////////////////////////////////////
// based on vertex-face collision pairs, build "collision tets"
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildVertexFaceCollisionTets(const VECTOR& velocity)
{
  // clear out old ones (not clear if this is smart or dumb)
  _vertexFaceCollisionTets.clear();
  _vertexFaceCollisionAreas.clear();

  // make a tet for each vertex-face pair
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const int vertexID = _vertexFaceCollisions[x].first;
    const int faceID = _vertexFaceCollisions[x].second;
    const VECTOR3I& face = _triangles[faceID];
    
    vector<VECTOR3> vs(4);
    for(int i = 0; i < 3; i++){
      vs[i] = _vertices[face[i]];
    }

    // build a tet with the correct vertex ordering
    VECTOR4I tet;
    tet[0] = vertexID;

    // cloth collisions: if face normal is aligned with collision
    // edge, invert face order. Otherwise feed face order normally

    VECTOR3 faceNorm = planeNormal(vs);
    VECTOR3 wTests = _vertices[vertexID] - _vertices[face[1]];
    if(wTests.dot(faceNorm) > 0)
    {
      tet[1] = face[2];
      tet[2] = face[1];
      tet[3] = face[0];
    }
    else
    {
      tet[1] = face[0];
      tet[2] = face[1];
      tet[3] = face[2];
    }
    // get the rest area of the triangle
    vector<VECTOR3> restFace(3);
    restFace[0] = _restVertices[face[0]]; 
    restFace[1] = _restVertices[face[1]];
    restFace[2] = _restVertices[face[2]];
    const REAL restFaceArea = triangleArea(restFace);

    const REAL restVertexArea = _restOneRingAreas[vertexID];

    // store everything
    _vertexFaceCollisionTets.push_back(tet);
    assert(restFaceArea >= 0.0);
    assert(restVertexArea >= 0.0);
    _vertexFaceCollisionAreas.push_back(restFaceArea + restVertexArea);
  }
}

///////////////////////////////////////////////////////////////////////
// compute collision forces using collision tets
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeVertexFaceCollisionForces() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++)
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    const VECTOR12 force = -_vertexFaceCollisionAreas[i] * _vertexFaceEnergy->gradient(vs);
    perElementForces[i] = force;

#if ENABLE_DEBUG_TRAPS
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " NaN in collision tet: " << i << endl;
      for (int j = 0; j < 4; j++)
        cout << " v" << j << ": " << vs[j].transpose() << endl;
      cout << " gradient: " << endl << _vertexFaceEnergy->gradient(vs) << endl;
    }
#endif
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const VECTOR12& tetForce = perElementForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}


///////////////////////////////////////////////////////////////////////
// compute edge-edge collision energy using x-based formulation
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeEdgeEdgeCollisionEnergy() const
{
  TIMER functionTimer(__FUNCTION__);

  REAL finalEnergy = 0.0;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;

    const REAL psi = _edgeEdgeEnergy->psi(vs,a,b);
    finalEnergy += _edgeEdgeCollisionAreas[i] * psi;
  }

  return finalEnergy;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision forces using x-based formulation
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeEdgeEdgeCollisionForces() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;

#if ADD_EDGE_EDGE_PENETRATION_BUG
    const VECTOR12 force = -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->gradient(vs,a,b);
#else
    bool intersect = _edgeEdgeIntersections[i];
    REAL area = _edgeEdgeCollisionAreas[i];
    const VECTOR12 force = (!intersect) ? -area * _edgeEdgeEnergy->gradient(vs,a,b)
                                                        : -area * _edgeEdgeEnergy->gradientNegated(vs,a,b);
#endif
    perElementForces[i] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];
    const VECTOR12& edgeForce = perElementForces[i];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * vertexIndices[x];
      assert((int)index < DOFs);
      forces[index]     += edgeForce[3 * x];
      forces[index + 1] += edgeForce[3 * x + 1];
      forces[index + 2] += edgeForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision Hessians using x-based formulation
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeEdgeEdgeCollisionClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;
#if ADD_EDGE_EDGE_PENETRATION_BUG
    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessian(vs,a,b);
#else
    const MATRIX12 H = (!_edgeEdgeIntersections[i]) ? -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessian(vs,a,b)
                                                    : -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessianNegated(vs,a,b);
#endif
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perElementHessians[i];
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int y = 0; y < 4; y++)
    {
      int yVertex = vertexIndices[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = vertexIndices[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }
 
  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the collision force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeVertexFaceCollisionClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++)
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    const MATRIX12 H = -_vertexFaceCollisionAreas[i] * _vertexFaceEnergy->clampedHessian(vs);
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

#if 0
///////////////////////////////////////////////////////////////////////
// find out how close all the edges are initially
//
// for each pair of _edges, what _collisionEps should we use? 
// If they started out closer than _collisionEps, then we need to set 
// a smaller tolerance.
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeEdgeRestDistance()
{
  TIMER functionTimer(__FUNCTION__);
  cout << " Computing edge-edge rest distances ... " << flush;
  _edgeEdgeRestDistance.clear();

  // hit all the pairs
  for (unsigned int y = 0; y < _edges.size(); y++)
  {
    // cache the vertex positions
    const VECTOR2I yEdge = _edges[y];
    const VECTOR3& v0 = _vertices[yEdge[0]];
    const VECTOR3& v1 = _vertices[yEdge[1]];

    // Wonder if meshes will get big enough that this needs to
    // be upgraded to 64-bit? Would be a good problem to have ...
    unsigned int yFlat = yEdge[0] + yEdge[1] * _edges.size();

    for (unsigned int x = 0; x < _edges.size(); x++)
    {
      // don't compare against yourself
      if (x == y) continue;

      const VECTOR2I xEdge = _edges[x];
      unsigned int xFlat = xEdge[0] + xEdge[1] * _edges.size();

      const VECTOR3& v2 = _vertices[xEdge[0]];
      const VECTOR3& v3 = _vertices[xEdge[1]];
    
      // do the actual distance test
      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);  
      const REAL distance = (innerPoint - outerPoint).norm();

      // store the pair
      pair<unsigned int,unsigned int> edgeEdge(xFlat, yFlat);

      // hash the distance
      _edgeEdgeRestDistance[edgeEdge] = distance;
    }
  }
  cout << "done. " << endl;
}
#endif

///////////////////////////////////////////////////////////////////////
// computes updated displacement from rest pose after processing CCDs
// Based off of Eberle's pseudocode description of Fizzt
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeCCDResponsePosition(const int maxIters)
{
  //First check if any CCDs are anywhere
  bool collisionFree = true;
  bool allRegionsResolved = false;
  if (_vertexFaceCCDs.size() + _edgeEdgeCCDs.size() != 0){
    collisionFree = false;
    gatherCCDs();
    // gatherCCDRegions();
  }
}

///////////////////////////////////////////////////////////////////////
// compute whether one vertex is inside the vertex one right of another
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceVertexOneRings()
{
  _insideSurfaceVertexOneRing.clear();
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    // gonna be lazy for a moment here
    const VECTOR2I edge = _edges[x];
    _insideSurfaceVertexOneRing[pair<int, int>(edge[0], edge[1])] = true;
    _insideSurfaceVertexOneRing[pair<int, int>(edge[1], edge[0])] = true;
  }
}

///////////////////////////////////////////////////////////////////////
// set collision eps to something new
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCollisionEps(const REAL& eps)
{
  _collisionEps = eps;
  _vertexFaceEnergy->eps() = eps;
  _edgeEdgeEnergy->setEps(eps);
}

///////////////////////////////////////////////////////////////////////
// set collision stiffness to something new
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCollisionStiffness(const REAL& mu)
{
  _vertexFaceEnergy->mu() = mu;
  _edgeEdgeEnergy->mu() = mu;
}

#if 0
///////////////////////////////////////////////////////////////////////
// compute the normal of the surface triangle at 
// _triangles[triangleID];
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::surfaceTriangleNormal(const int triangleID) const
{
  assert(triangleID < (int)_triangles.size());

  const VECTOR3I& vertexIDs = _triangles[triangleID];
  const VECTOR3& v0 = _vertices[vertexIDs[0]];
  const VECTOR3& v1 = _vertices[vertexIDs[1]];
  const VECTOR3& v2 = _vertices[vertexIDs[2]];

  const VECTOR3& e0 = v1 - v0;
  const VECTOR3& e1 = v2 - v0;

  return e0.cross(e1).normalized();
}

///////////////////////////////////////////////////////////////////////
// set collision pairs (for replays)
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                                 const vector<pair<int,int> >& edgeEdge)
{
  _vertexFaceCollisions = vertexFace;
  _edgeEdgeCollisions = edgeEdge;
}

///////////////////////////////////////////////////////////////////////
// compute the dihedral angle between two surface faces
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::surfaceFaceDihedralAngle(const int surfaceID0, const int surfaceID1) const
{
  const VECTOR4I tet = buildSurfaceFlap(surfaceID0, surfaceID1);

  // let's do some cross products ...
  ///////////////////////////
  //                       //                                                                         
  //           1           //
  //                       //
  //           o           //
  //          /|\          //
  //         / | \         //
  //        /  |  \        //
  //    0  o   |   o  3    //
  //        \  |  /        //
  //         \ | /         //
  //          \|/          //
  //           o           //
  //                       //
  //           2           //
  //                       //
  ///////////////////////////
  const VECTOR3& v0 = _vertices[tet[0]];
  const VECTOR3& v1 = _vertices[tet[1]];
  const VECTOR3& v2 = _vertices[tet[2]];
  const VECTOR3& v3 = _vertices[tet[3]];

  const VECTOR3 e20 = v2 - v0;
  const VECTOR3 e10 = v1 - v0;
  const VECTOR3 n0 = e20.cross(e10) / (e20 - e10).norm();

  const VECTOR3 e13 = v1 - v3;
  const VECTOR3 e23 = v2 - v3;
  const VECTOR3 n1 = e13.cross(e23) / (e13 - e23).norm();

  const VECTOR3 e12 = (v1 - v2) / (v1 - v2).norm();

  const REAL sinTheta = (n0.cross(n1)).dot(e12);
  const REAL cosTheta = n0.dot(n1);

  return atan2(sinTheta, cosTheta);
}
#endif

///////////////////////////////////////////////////////////////////////
// build a consistent tet/flap ordering from two surface triangles
///////////////////////////////////////////////////////////////////////
VECTOR4I TRIANGLE_MESH::buildFlap(const int triangle0, const int triangle1) const
{
  // they're legal indices, right?
  assert(triangle0 >= 0);
  assert(triangle1 >= 0);
  assert(triangle0 < (int)_triangles.size());
  assert(triangle1 < (int)_triangles.size());

  // they are neighbors, right?
  assert(areTriangleNeighbors(triangle0, triangle1));

  const VECTOR3I f0 = _triangles[triangle0];
  const VECTOR3I f1 = _triangles[triangle1];

  int firstMatch = -1;
  int secondMatch = -1;
  int unmatched0 = -1;
  int unmatched1 = -1;

  // find the tet indices for the first face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i0 = f0[x];

    bool matchFound = false;
    for (int y = 0; y < 3; y++)
    {
      // see if it matches
      if (i0 == f1[y])
      {
        // which match is it?
        if (firstMatch == -1)
          firstMatch = i0;
        else
          secondMatch = i0;
        
        matchFound = true;
      }
    }

    if (!matchFound)
      unmatched0 = i0;
  }

  // find the unmatched vertex from the second face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i1 = f1[x];
    if (i1 != firstMatch && i1 != secondMatch)
      unmatched1 = i1;
  }

  // we did find one, right?
  assert(unmatched1 != -1);

  // build a tet/flap
  ///////////////////////////
  //                       //                                                                         
  //           1           //
  //                       //
  //           o           //
  //          /|\          //
  //         / | \         //
  //        /  |  \        //
  //    0  o   |   o  3    //
  //        \  |  /        //
  //         \ | /         //
  //          \|/          //
  //           o           //
  //                       //
  //           2           //
  //                       //
  ///////////////////////////
  VECTOR4I tet;
  tet[0] = unmatched0;
  tet[1] = secondMatch;
  tet[2] = firstMatch;
  tet[3] = unmatched1;

  return tet;
}

///////////////////////////////////////////////////////////////////////
// compute flaps and their areas
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFlaps()
{
  _flaps.clear();
  for (unsigned int x = 0; x < _edgeTriangleNeighbors.size(); x++)
  {
    const VECTOR2I& edgeTriangles = _edgeTriangleNeighbors[x];
    if (edgeTriangles[0] == -1 || edgeTriangles[1] == -1) continue;

    const VECTOR4I flap = buildFlap(edgeTriangles[0], 
                                    edgeTriangles[1]);
    _flaps.push_back(flap);

    const REAL area = _restTriangleAreas[edgeTriangles[0]] + 
                      _restTriangleAreas[edgeTriangles[1]];
    _restFlapAreas.push_back(area);
  }
}

///////////////////////////////////////////////////////////////////////
// get mass-weighted global translation
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::getTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL areaSum = 0;
  assert(_vertices.size() == _restOneRingAreas.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    areaSum += _restOneRingAreas[x];
    vertexSum += _vertices[x] * _restOneRingAreas[x];
  }

  return vertexSum * (1.0 / areaSum);
}

///////////////////////////////////////////////////////////////////////
// get volume-weighted global translation, for the rest state
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::getRestTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL areaSum = 0;
  assert(_vertices.size() == _restOneRingAreas.size());
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    areaSum += _restOneRingAreas[x];
    vertexSum += _restVertices[x] * _restOneRingAreas[x];
  }

  return vertexSum * (1.0 / areaSum);
}

///////////////////////////////////////////////////////////////////////
// are these two surface triangles neighbors?
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::areTriangleNeighbors(const int id0, const int id1) const
{
  assert(_triangleNeighbors.size() > 0);
  assert(id0 < (int)_triangleNeighbors.size());
  assert(id1 < (int)_triangleNeighbors.size());

  const VECTOR3I neighbors0 = _triangleNeighbors[id0];

  for (int x = 0; x < 3; x++)
    if (neighbors0[x] == id1)
      return true;

  return false;
}

#if 0
///////////////////////////////////////////////////////////////////////
// get Procrustes-style global rotation using Eqn. 7 and surrounding
// text from Muller et al's "Meshless Deformations Based on Shape 
// Matching" from SIGGRAPH 2005
///////////////////////////////////////////////////////////////////////
MATRIX3 TRIANGLE_MESH::getRotation() const
{
  // trying to follow Muller's notation here
  const VECTOR3 x_cm0 = getRestTranslation();
  const VECTOR3 x_cm  = getTranslation();

  // left matrix in Eqn. 7 of Muller's paper
  MATRIX3 Apq;
  Apq.setZero();
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    const VECTOR3 p = _vertices[x] - x_cm;
    const VECTOR3 q = _restVertices[x] - x_cm0;
    Apq += _restOneRingVolumes[x] * (p * q.transpose());
  }

  // get the rotation
  MATRIX3 R,S;
  polarDecomposition(Apq, R, S);
  return R;
}
#endif

///////////////////////////////////////////////////////////////////////
// see if a current surface triangle has been crushed to degeneracy
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::surfaceTriangleIsDegenerate(const int surfaceTriangleID)
{
  assert(surfaceTriangleID >= 0);
  assert(surfaceTriangleID < (int)_triangles.size());

  // get the rest area
  vector<VECTOR3> vertices(3);
  vertices[0] = _restVertices[_triangles[surfaceTriangleID][0]];
  vertices[1] = _restVertices[_triangles[surfaceTriangleID][1]];
  vertices[2] = _restVertices[_triangles[surfaceTriangleID][2]];
  const REAL restArea = triangleArea(vertices);

  // get the deformed area
  vertices[0] = _vertices[_triangles[surfaceTriangleID][0]];
  vertices[1] = _vertices[_triangles[surfaceTriangleID][1]];
  vertices[2] = _vertices[_triangles[surfaceTriangleID][2]];
  const REAL deformedArea = triangleArea(vertices);

  const REAL relativeArea = deformedArea / restArea;

  const REAL degeneracyEps = 1e-4;
  if (relativeArea < degeneracyEps) return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute which vertices are attached to inverted tets
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeInvertedVertices()
{
  TIMER functionTimer(__FUNCTION__);
  // first set them all to false
  _invertedVertices.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _invertedVertices[x] = false;

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // if the tet is not inverted, move on
    if (invariant3(_Fs[x]) > 0.0)
      continue;

    // if tet is inverted, tags all its vertices
    for (int y = 0; y < 3; y++)
      _invertedVertices[_triangles[x][y]] = true;
  }

  int totalInverted = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (_invertedVertices[x])
      totalInverted++;

  //cout << " Total inverted vertices: " << totalInverted << endl;
}

}
