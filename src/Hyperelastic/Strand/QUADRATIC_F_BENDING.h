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
// July 15, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STRAND_QUADRATIC_F_BENDING_H
#define STRAND_QUADRATIC_F_BENDING_H

#include "TAN_BENDING.h"


namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
class QUADRATIC_F_BENDING : public TAN_BENDING
{
public:
  QUADRATIC_F_BENDING(const REAL& mu);
  virtual ~QUADRATIC_F_BENDING() { };

  // Computes the strain energy density
  // E is the matrix of edge vectors, stacked column-wise
  virtual REAL psi(const MATRIX3x2& E, const REAL& theta0) const override;
  virtual MATRIX3x2 PK1(const MATRIX3x2& E, const REAL& theta0) const override;
  virtual MATRIX6 hessian(const MATRIX3x2& E, const REAL& theta0) const override;
  virtual MATRIX6 clampedHessian(const MATRIX3x2& E, const REAL& theta0) const override;

  virtual REAL psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const override;
  virtual MATRIX3x2 PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const override;
  virtual MATRIX6 hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const override;
  virtual MATRIX6 clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const override;

  // The name of the material
  virtual std::string name() const override;

protected:
  // gradient computation, once and for all, with divide-by-zero guard
  MATRIX3x2 gradient(const MATRIX3x2& E) const;
};

}
}

#endif
