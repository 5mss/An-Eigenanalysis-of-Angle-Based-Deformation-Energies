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
// This is a file in the ANGLE library
//
// May 26, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHELL_QUADRATIC_F_BENDING_H
#define SHELL_QUADRATIC_F_BENDING_H

#include "BENDING.h"

#define USING_UNFILTERED 0
#define USING_GAUSS_NEWTON_BENDING 0
#define USING_BRUTE_FORCE_CLAMP 0
#define USING_RANK4_CORRECTION 1

namespace HOBAK {
namespace SHELL {

class QUADRATIC_F_BENDING : public BENDING
{
public:
  QUADRATIC_F_BENDING(const REAL& mu);
  virtual ~QUADRATIC_F_BENDING() { };

  // Computes the strain energy density
  virtual REAL psi(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // position-based gradient
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // position-based Hessian
  virtual MATRIX12 hessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // filtered hessian
  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // The name of the material
  virtual std::string name() const;

  // Compute the bending angle based on a rest flap.
  virtual REAL restAngle(const std::vector<VECTOR3>& restFlap) const;

protected:

};

}
}

#endif