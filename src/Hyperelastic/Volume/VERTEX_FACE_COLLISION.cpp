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
#include "VERTEX_FACE_COLLISION.h"
#include "MATRIX_UTIL.h"
#include "COLLISION_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

VERTEX_FACE_COLLISION::VERTEX_FACE_COLLISION(const REAL& mu, const REAL& eps) :
    _mu(mu), _eps(eps)
{
}

std::string VERTEX_FACE_COLLISION::name() const
{
  return "VERTEX_FACE_COLLISION";
}

///////////////////////////////////////////////////////////////////////
// convert the 12-vector in a way that imposes a consistent tet 
// ordering for vertices and edges
//
// assumes v0, x[0,1,2] is the collision vertex
///////////////////////////////////////////////////////////////////////
void VERTEX_FACE_COLLISION::getVerticesAndEdges(const VECTOR12& x,
                                                vector<VECTOR3>& v,
                                                vector<VECTOR3>& e)
{
  v.resize(4);
  for (int i = 0; i < 4; i++)
  {
    v[i][0] = x[i * 3];
    v[i][1] = x[i * 3 + 1];
    v[i][2] = x[i * 3 + 2];
  }

  e.resize(3);
  e[0] = v[3] - v[2];
  e[1] = v[0] - v[2];
  e[2] = v[1] - v[2];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL VERTEX_FACE_COLLISION::psi(const vector<VECTOR3>& vertices) const
{
  return psi(flattenVertices(vertices));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 VERTEX_FACE_COLLISION::gradient(const vector<VECTOR3>& vertices) const
{
  return gradient(flattenVertices(vertices));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_COLLISION::hessian(const vector<VECTOR3>& vertices) const
{
  return hessian(flattenVertices(vertices));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL VERTEX_FACE_COLLISION::psi(const VECTOR12& x) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);

  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // get the spring length, non-zero rest-length
  const VECTOR3 tvf = v[0] - v[2];
  REAL springLength = (n.dot(tvf) - _eps);
  return _mu * springLength * springLength;
}

///////////////////////////////////////////////////////////////////////
// gradient of spring length, n' * (v[0] - v[2])
///////////////////////////////////////////////////////////////////////
VECTOR12 VERTEX_FACE_COLLISION::springLengthGradient(const vector<VECTOR3>& v,
                                                     const vector<VECTOR3>& e,
                                                     const VECTOR3& n)
{
  const MATRIX3x12 nPartial = normalGradientVF(e);
  const VECTOR3 tvf = v[0] - v[2];

  MATRIX3x12 tvfPartial;
  tvfPartial.setZero();
  tvfPartial(0,0) = tvfPartial(1,1) = tvfPartial(2,2) = 1.0;
  tvfPartial(0,6) = tvfPartial(1,7) = tvfPartial(2,8) = -1.0;

  //f = nPartial' * (v2 - v0) + tvfPartial' * n;
  return nPartial.transpose() * tvf + tvfPartial.transpose() * n;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 VERTEX_FACE_COLLISION::gradient(const VECTOR12& x) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();
  
  // get the spring length, non-zero rest-length
  const VECTOR3 tvf = v[0] - v[2];
  REAL springLength = n.dot(tvf) - _eps;
  return 2.0 * _mu * springLength * springLengthGradient(v,e,n);
}

///////////////////////////////////////////////////////////////////////
// hessian of spring length, n' * (v[0] - v[2])
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_COLLISION::springLengthHessian(const vector<VECTOR3>& v,
                                                    const vector<VECTOR3>& e,
                                                    const VECTOR3& n)
{
  const VECTOR3 tvf = v[0] - v[2];

  MATRIX3x12 tvfPartial;
  tvfPartial.setZero();
  tvfPartial(0,0) = tvfPartial(1,1) = tvfPartial(2,2) = 1.0;
  tvfPartial(0,6) = tvfPartial(1,7) = tvfPartial(2,8) = -1.0;

  //% mode-3 contraction
  //[nx ny nz] = normal_hessian(x);
  //final = nx * tvf(1) + ny * tvf(2) + nz * tvf(3);
  const vector<MATRIX12> normalH = normalHessianVF(e);
  const MATRIX12 contracted = tvf[0] * normalH[0] + tvf[1] * normalH[1] + 
                              tvf[2] * normalH[2];
  
  const MATRIX3x12 nGrad = normalGradientVF(e);

  //product = nGrad' * vGrad;
  const MATRIX12 product = nGrad.transpose() * tvfPartial;

  return contracted + product + product.transpose();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_COLLISION::hessian(const VECTOR12& x) const
{
  // convert to vertices and edges
  vector<VECTOR3> v;
  vector<VECTOR3> e;
  getVerticesAndEdges(x, v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // get the spring length, non-zero rest-length
  const VECTOR3 tvf = v[0] - v[2];
  const REAL springLength = n.dot(tvf) - _eps;

  // ndotGrad    = ndot_gradient(x);
  const VECTOR12 gvf = springLengthGradient(v,e,n);

  // ndotHessian = ndot_hessian(x);
  const MATRIX12 springLengthH = springLengthHessian(v,e,n);
  
  // final = 2 * k * (ndotGrad * ndotGrad' + ndot * ndotHessian);
  return 2.0 * _mu * (gvf * gvf.transpose() + 
                      springLength * springLengthH);
  // Gauss-Newton approximation
  //return 2.0 * _mu * (gvf * gvf.transpose()); 
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_COLLISION::clampedHessian(const VECTOR12& x) const
{
  return clampEigenvalues(hessian(x));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 VERTEX_FACE_COLLISION::clampedHessian(const vector<VECTOR3>& v) const
{
  return clampedHessian(flattenVertices(v));
}

} // VOLUME
} // HOBAK
