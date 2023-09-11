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
#ifndef VOLUME_EDGE_COLLISION_H
#define VOLUME_EDGE_COLLISION_H

#include "SETTINGS.h"
#include <vector>

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is (sort of) the edge-edge collision energy described in a variety of
// places:
//
// "Dynamic deformables: implementation and production practicalities"
// Kim and Eberle, 2020, Chapter 11.4
//
// "Robust treatement of simultaneous contacts"
// Harmon et al., 2008, Section 2
//
// "Collision and self-collision handling in cloth model dedicated to design garments"
// Provot 1997, Section 3.3
//
// I haven't actually seen it written down as an energy anywhere. The Harmon
// paper comes closest by posing an unsquared constraint.
//
// The energy is extremely similar to the vertex-face collision energy from:
//
// "Efficient elasticity for character skinning with contact and collisions"
// McAdams et al. 2011.
//
// The normal gradient and Hessian stuff works out to almost the same thing,
// up to a few constants
//
// vertex and coordinate ordering is:
//
//    vector<VECTOR3> vs(4);
//    vs[0] = _vertices[edge0[0]];
//    vs[1] = _vertices[edge0[1]];
//    vs[2] = _vertices[edge1[0]];
//    vs[3] = _vertices[edge1[1]];
//    
//    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
//    const VECTOR2& b = _edgeEdgeCoordinates[i].second;
//
///////////////////////////////////////////////////////////////////////////////////
class EDGE_COLLISION
{
public:
  EDGE_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  virtual ~EDGE_COLLISION() {};

  // get the strain energy
  REAL psi(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual REAL psi(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;
  
  // negated version of strain energy, if it's already intersecting
  REAL psiNegated(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual REAL psiNegated(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  VECTOR12 gradient(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual VECTOR12 gradient(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;
  
  // negated version of gradient, if it's already intersecting
  VECTOR12 gradientNegated(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual VECTOR12 gradientNegated(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;

  MATRIX12 hessian(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual MATRIX12 hessian(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;
  
  MATRIX12 hessianNegated(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual MATRIX12 hessianNegated(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;

  MATRIX12 clampedHessian(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual MATRIX12 clampedHessian(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;
  
  MATRIX12 clampedHessianNegated(const std::vector<VECTOR3>& v, const VECTOR2& a, const VECTOR2& b) const;
  virtual MATRIX12 clampedHessianNegated(const VECTOR12& x, const VECTOR2& a, const VECTOR2& b) const;
  
  virtual std::string name() const;

  const REAL& mu() const  { return _mu; };
  const REAL& eps() const { return _eps; };
  REAL& mu()  { return _mu; };
  //REAL& eps() { return _eps; };
  virtual void setEps(const REAL& eps) { _eps = eps; };

  // convert vertices to edges before determining if the two
  // edges are nearly parallel. this is needed by unit tests to
  // determine if we should give up on the current test.
  static bool nearlyParallelVertices(const std::vector<VECTOR3>& x);

protected:
  // convert the 12-vector in a way that imposes a consistent tet 
  // ordering for vertices and edges
  static void getVerticesAndEdges(const VECTOR12& x, 
                                  std::vector<VECTOR3>& v, 
                                  std::vector<VECTOR3>& e);

  // gradient of spring length, n' * (va - vb)
  static VECTOR12 springLengthGradient(const std::vector<VECTOR3>& e,
                                       const VECTOR3& n,
                                       const VECTOR3& diff,
                                       const VECTOR2& a,
                                       const VECTOR2& b);

  // hessian of spring length, n' * (va - vb)
  static MATRIX12 springLengthHessian(const std::vector<VECTOR3>& e,
                                      const VECTOR3& n,
                                      const VECTOR3& diff,
                                      const VECTOR2& a,
                                      const VECTOR2& b);

  // partial of (va - vb)
  static MATRIX3x12 vDiffPartial(const VECTOR2& a, const VECTOR2& b);

  // are the two edges nearly parallel?
  static bool nearlyParallelEdges(const std::vector<VECTOR3>& e);

  // collision stiffness
  REAL _mu;

  // collision epsilon -- how far apart should we push things?
  REAL _eps;
};

} // VOLUME
} // HOBAK

#endif
