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
#ifndef VOLUME_VERTEX_FACE_COLLISION_H
#define VOLUME_VERTEX_FACE_COLLISION_H

#include "SETTINGS.h"
#include <vector>

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is (sort of) the collision energy described in the SIGGRAPH 2011 paper 
//
// "Efficient elasticity for character skinning with contact and collisions"
//
// by Alex McAdams, Yongning Zhu, Andrew Selle, Mark Empey, Rasmus Tamstorf,
// Joseph Teran, and Eftychios Sifakis
//
// This energy is from the last paragraph of Section 7, but the x_s term there
// appears to be fixed; it's a proxy point sprinkled onto the body surface in
// the style of "Hybrid simulation of deformable solids" by Sifakis et al. 2007.
//
// The formulation implemented here is the more general case where the 
// collision point on the triangle slides along with the collision vertex. As 
// such, the second \alpha \II term in their energy disappears. I have a document 
// about this somewhere, if you want more details. (T. Kim 8/19/21)
//
// I added it to the course notes! It's now in the collision chapter, under
// the subsection "Another (Better? Identical?) Vertex-Face Energy"
// (T. Kim 1/7/22)
///////////////////////////////////////////////////////////////////////////////////
class VERTEX_FACE_COLLISION
{
public:
  VERTEX_FACE_COLLISION(const REAL& mu, const REAL& eps = 0.0);
  virtual ~VERTEX_FACE_COLLISION() {};

  // get the strain energy
  virtual REAL psi(const std::vector<VECTOR3>& v) const;
  virtual REAL psi(const VECTOR12& x) const;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& v) const;
  virtual VECTOR12 gradient(const VECTOR12& x) const;

  virtual std::string name() const;

  virtual MATRIX12 hessian(const std::vector<VECTOR3>& v) const;
  virtual MATRIX12 hessian(const VECTOR12& x) const;

  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& v) const;
  virtual MATRIX12 clampedHessian(const VECTOR12& x) const;

  const REAL& mu() const  { return _mu; };
  const REAL& eps() const { return _eps; };
  REAL& mu()  { return _mu; };
  REAL& eps() { return _eps; };
    
protected:
  // convert the 12-vector in a way that imposes a consistent tet 
  // ordering for vertices and edges
  static void getVerticesAndEdges(const VECTOR12& x, 
                                  std::vector<VECTOR3>& v, 
                                  std::vector<VECTOR3>& e);

  // gradient of spring length, n' * (v[0] - v[2])
  static VECTOR12 springLengthGradient(const std::vector<VECTOR3>& v,
                                       const std::vector<VECTOR3>& e,
                                       const VECTOR3& n);

  // hessian of spring length, n' * (v[0] - v[2])
  static MATRIX12 springLengthHessian(const std::vector<VECTOR3>& v,
                                      const std::vector<VECTOR3>& e,
                                      const VECTOR3& n);

  // collision stiffness
  REAL _mu;

  // collision epsilon -- how far apart should we push things?
  REAL _eps;
};

} // VOLUME
} // ANGLE

#endif
