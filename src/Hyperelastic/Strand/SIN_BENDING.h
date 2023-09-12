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
// November 8, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STRAND_SIN_BENDING_H
#define STRAND_SIN_BENDING_H

#include "SETTINGS.h"

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
class SIN_BENDING
{
public:
  SIN_BENDING();
  virtual ~SIN_BENDING() { };

  // The name of the material
  virtual std::string name() const;

  REAL psi(const MATRIX2& B, const VECTOR2& kappa, const VECTOR2& kappaBar) const;

  // ugh, so verbose.
  MATRIX11 hessian(const MATRIX3x2& E, 
                   const MATRIX2& B, 
                   const MATRIX3x2& M0, 
                   const MATRIX3x2& M1, 
                   const MATRIX3x2& M, 
                   const VECTOR2& kappa, 
                   const VECTOR2& kappaBar) const;
  VECTOR11 gradient(const MATRIX3x2& E, 
                    const MATRIX2& B, 
                    const MATRIX3x2& M0, 
                    const MATRIX3x2& M1, 
                    const MATRIX3x2& M, 
                    const VECTOR2& kappa, 
                    const VECTOR2& kappaBar) const;
  virtual MATRIX3x2 edgeGradient(const MATRIX3x2& E, 
                                 const MATRIX2& B, 
                                 const MATRIX3x2& M, 
                                 const VECTOR2& kappa, 
                                 const VECTOR2& kappaBar) const;
  VECTOR2 twistGradient(const MATRIX3x2& E, 
                        const MATRIX2& B, 
                        const MATRIX3x2& M0, 
                        const MATRIX3x2& M1, 
                        const VECTOR2& kappa, 
                        const VECTOR2& kappaBar) const;
  virtual VECTOR3 binormal(const VECTOR3& e0, const VECTOR3& e1) const;
  VECTOR3 binormal(const MATRIX3x2& E) const;

  MATRIX6 edgeHessian(const MATRIX3x2& E,
                      const MATRIX2& B,
                      const MATRIX3x2& M,
                      const VECTOR2& kappaBar) const;
  MATRIX2 twistHessian(const MATRIX3x2& E, 
                       const MATRIX2& B, 
                       const MATRIX3x2& M0, 
                       const MATRIX3x2& M1, 
                       const VECTOR2& kappa, 
                       const VECTOR2& kappaBar) const;
  std::pair<MATRIX3x2, MATRIX3x2> mixedHessian(const MATRIX3x2& E, 
                                               const MATRIX2& B, 
                                               const MATRIX3x2& M0, 
                                               const MATRIX3x2& M1, 
                                               const MATRIX3x2& M, 
                                               const VECTOR2& kappa, 
                                               const VECTOR2& kappaBar) const;

  virtual std::pair<MATRIX3,MATRIX3> binormalGradient(const MATRIX3x2& E) const;
private:
  virtual std::tuple<MATRIX3,MATRIX3,MATRIX3> binormalHessian00(const MATRIX3x2& E) const;
  virtual std::tuple<MATRIX3,MATRIX3,MATRIX3> binormalHessian11(const MATRIX3x2& E) const;
  virtual std::tuple<MATRIX3,MATRIX3,MATRIX3> binormalHessian01(const MATRIX3x2& E) const;
  static MATRIX6x9 dedx();
};

}
}

#endif
