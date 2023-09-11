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
///////////////////////////////////////////////////////////////////////////////////////////////////
// This is a file in the HOBAK library
//
// July 11, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STRAND_QUADRATIC_STRETCHING_H
#define STRAND_QUADRATIC_STRETCHING_H

#include "STRETCHING.h"

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////
// the strand geometry is assumed to be:
//
//      p0                          p1
//       o--------------------------o
//
// the deformation gradient "f" is assumed to be:
//
//  f = p1 - p0
//      -------
//     |p1 - p0|    (the magnitude / 2-norm of the difference)
//
// this roughly corresponds to the other treatments we have seen of F = D_s * D_m^{-1},
// where we treat:
//
//    D_s = p_1 - p_0
//    D_m = |p_1 - p_0|
//
///////////////////////////////////////////////////////////////////////////////////////////////////
class QUADRATIC_STRETCHING : public STRETCHING
{
public:
  QUADRATIC_STRETCHING(const REAL& mu);
  virtual ~QUADRATIC_STRETCHING() { };

  // Computes the strain energy density
  virtual REAL psi(const VECTOR3& f) const override;

  // f-based gradient
  virtual VECTOR3 PK1(const VECTOR3& f) const override;

  // f-based Hessian
  virtual MATRIX3 hessian(const VECTOR3& f) const override;
 
  // f-based clamped Hessian
  virtual MATRIX3 clampedHessian(const VECTOR3& f) const override;

  // position-based gradient
  virtual VECTOR6 spatialGradient(const std::vector<VECTOR3>& p, const REAL& dmInv) const override;

  // position-based Hessian
  virtual MATRIX6 spatialHessian(const std::vector<VECTOR3>& p, const REAL& dmInv) const override;

  // position-based clamped Hessian
  virtual MATRIX6 spatialClampedHessian(const std::vector<VECTOR3>& p, const REAL& dmInv) const override;

  // The name of the material
  virtual std::string name() const override;
};

}
}

#endif
