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

#include "SIN_BENDING.h"
#include "MATRIX_UTIL.h"

#include <iostream>

using namespace std;

namespace HOBAK {
namespace STRAND {

SIN_BENDING::SIN_BENDING()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string SIN_BENDING::name() const
{ 
  return std::string("SIN_BENDING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL SIN_BENDING::psi(const MATRIX2& B, const VECTOR2& kappa, const VECTOR2& kappaBar) const
{
  //return (kappa - kappaBar).transpose() * B * (kappa - kappaBar);
  return 0.5 * (kappa - kappaBar).transpose() * B * (kappa - kappaBar);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 SIN_BENDING::edgeGradient(const MATRIX3x2& E, 
                                    const MATRIX2& B, 
                                    const MATRIX3x2& M, 
                                    const VECTOR2& kappa,
                                    const VECTOR2& kappaBar) const
{
  // get the binormal gradient
  const std::pair<MATRIX3, MATRIX3> sinGradient = binormalGradient(E);
    
  // get the edge gradient
  const VECTOR3 rhs = M * B * (kappa - kappaBar);
  MATRIX3x2 gradient;
  //gradient.col(0) = 2.0 * sinGradient.first.transpose()  * rhs;
  //gradient.col(1) = 2.0 * sinGradient.second.transpose() * rhs;
  gradient.col(0) = sinGradient.first.transpose()  * rhs;
  gradient.col(1) = sinGradient.second.transpose() * rhs;

  return gradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR2 SIN_BENDING::twistGradient(const MATRIX3x2& E, 
                                   const MATRIX2& B, 
                                   const MATRIX3x2& M0, 
                                   const MATRIX3x2& M1, 
                                   const VECTOR2& kappa,
                                   const VECTOR2& kappaBar) const
{
  const VECTOR3 kb = binormal(E.col(0), E.col(1));
  VECTOR2 gradient;
  gradient[0] = -kb.transpose() * M0 * B * (kappa - kappaBar);
  gradient[1] = -kb.transpose() * M1 * B * (kappa - kappaBar);

  //return gradient;
  return 0.5 * gradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 SIN_BENDING::binormal(const MATRIX3x2& E) const
{
  return binormal(E.col(0), E.col(1));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 SIN_BENDING::binormal(const VECTOR3& e0, const VECTOR3& e1) const
{
  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();

  return e0.cross(e1) / (e0norm * e1norm);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<MATRIX3,MATRIX3> SIN_BENDING::binormalGradient(const MATRIX3x2& E) const
{
  const VECTOR3& e0 = E.col(0);
  const VECTOR3& e1 = E.col(1);
  const VECTOR3 kb = binormal(E);
  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();

  const MATRIX3 cross0 = crossProduct(e0);
  const MATRIX3 cross1 = crossProduct(e1);

  std::pair<MATRIX3, MATRIX3> G;

  const REAL coeff = 1.0 / (e0norm * e1norm);
  G.first  = -(1.0 / (e0norm * e0norm)) * kb * e0.transpose() - coeff * cross1;
  G.second = -(1.0 / (e1norm * e1norm)) * kb * e1.transpose() + coeff * cross0;
  return G;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple<MATRIX3,MATRIX3,MATRIX3> SIN_BENDING::binormalHessian00(const MATRIX3x2& E) const
{
  const VECTOR3& e0 = E.col(0);
  const VECTOR3& e1 = E.col(1);
  const VECTOR3 kb = binormal(E);
  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();

  const MATRIX3 cross1 = crossProduct(e1);
  const std::pair<MATRIX3,MATRIX3> G = binormalGradient(E);
  const MATRIX3& G0 = G.first;

  const VECTOR3 ex(1,0,0);
  const VECTOR3 ey(0,1,0);
  const VECTOR3 ez(0,0,1);

  const REAL denom = e0norm * e0norm;

  const MATRIX3 X = (1.0 / denom) * ((2.0 * e0[0] / denom) * (kb * e0.transpose()) - G0.col(0) * e0.transpose() - kb * ex.transpose() + e0[0] / (e0norm * e1norm) * cross1); 
  const MATRIX3 Y = (1.0 / denom) * ((2.0 * e0[1] / denom) * (kb * e0.transpose()) - G0.col(1) * e0.transpose() - kb * ey.transpose() + e0[1] / (e0norm * e1norm) * cross1); 
  const MATRIX3 Z = (1.0 / denom) * ((2.0 * e0[2] / denom) * (kb * e0.transpose()) - G0.col(2) * e0.transpose() - kb * ez.transpose() + e0[2] / (e0norm * e1norm) * cross1); 
  return std::make_tuple(X,Y,Z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple<MATRIX3,MATRIX3,MATRIX3> SIN_BENDING::binormalHessian11(const MATRIX3x2& E) const
{
  const VECTOR3& e0 = E.col(0);
  const VECTOR3& e1 = E.col(1);
  const VECTOR3 kb = binormal(E);
  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();

  const MATRIX3 cross0 = crossProduct(e0);
  const std::pair<MATRIX3,MATRIX3> G = binormalGradient(E);
  const MATRIX3& G1 = G.second;

  const VECTOR3 ex(1,0,0);
  const VECTOR3 ey(0,1,0);
  const VECTOR3 ez(0,0,1);

  const REAL denom = e1norm * e1norm;

  const MATRIX3 X = (1.0 / denom) * ((2.0 * e1[0] / denom) * (kb * e1.transpose()) - G1.col(0) * e1.transpose() - kb * ex.transpose() - e1[0] / (e0norm * e1norm) * cross0); 
  const MATRIX3 Y = (1.0 / denom) * ((2.0 * e1[1] / denom) * (kb * e1.transpose()) - G1.col(1) * e1.transpose() - kb * ey.transpose() - e1[1] / (e0norm * e1norm) * cross0); 
  const MATRIX3 Z = (1.0 / denom) * ((2.0 * e1[2] / denom) * (kb * e1.transpose()) - G1.col(2) * e1.transpose() - kb * ez.transpose() - e1[2] / (e0norm * e1norm) * cross0); 
  return std::make_tuple(X,Y,Z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::tuple<MATRIX3,MATRIX3,MATRIX3> SIN_BENDING::binormalHessian01(const MATRIX3x2& E) const
{
  const VECTOR3& e0 = E.col(0);
  const VECTOR3& e1 = E.col(1);
  const REAL e0norm = e0.norm();
  const REAL e1norm = e1.norm();

  const MATRIX3 cross1 = crossProduct(e1);
  const VECTOR3 ex(1,0,0);
  const VECTOR3 ey(0,1,0);
  const VECTOR3 ez(0,0,1);
  const MATRIX3 crossX = crossProduct(ex);
  const MATRIX3 crossY = crossProduct(ey);
  const MATRIX3 crossZ = crossProduct(ez);

  const std::pair<MATRIX3,MATRIX3> G = binormalGradient(E);
  const MATRIX3& G1 = G.second;

  const REAL denom = e1norm * e0norm;
  const REAL e0dot = e0norm * e0norm;
  const REAL e1dot = e1norm * e1norm;

  const MATRIX3 X = -(1.0 / e0dot) * G1.col(0) * e0.transpose() + e1[0] / (denom * e1dot) * cross1 - 1.0 / denom * crossX;
  const MATRIX3 Y = -(1.0 / e0dot) * G1.col(1) * e0.transpose() + e1[1] / (denom * e1dot) * cross1 - 1.0 / denom * crossY;
  const MATRIX3 Z = -(1.0 / e0dot) * G1.col(2) * e0.transpose() + e1[2] / (denom * e1dot) * cross1 - 1.0 / denom * crossZ;
  return std::make_tuple(X,Y,Z);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6x9 SIN_BENDING::dedx()
{
  MATRIX6x9 result;
  result.setZero();
  result.block<3,3>(0,0) = -MATRIX3::Identity();
  result.block<3,3>(0,3) =  MATRIX3::Identity();
  result.block<3,3>(3,3) = -MATRIX3::Identity();
  result.block<3,3>(3,6) =  MATRIX3::Identity();
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 SIN_BENDING::edgeHessian(const MATRIX3x2& E,
                                 const MATRIX2& B,
                                 const MATRIX3x2& M,
                                 const VECTOR2& kappaBar) const
{
  const VECTOR3 ks = binormal(E);
  const VECTOR2 kappa = M.transpose() * ks;
  const VECTOR3 rhs = M * B * (kappa - kappaBar);
  const pair<MATRIX3,MATRIX3> G = binormalGradient(E);
  const MATRIX3 E0 = G.first;
  const MATRIX3 E1 = G.second;
 
  MATRIX3 term1, term2;

  const std::tuple<MATRIX3, MATRIX3, MATRIX3> BH00 = binormalHessian00(E);
  term1.col(0) = std::get<0>(BH00).transpose() * rhs;
  term1.col(1) = std::get<1>(BH00).transpose() * rhs;
  term1.col(2) = std::get<2>(BH00).transpose() * rhs;
  //term1 *= 2.0;
  //term2 = 2.0 * E0.transpose() * M * B * M.transpose() * E0;
  term2 = E0.transpose() * M * B * M.transpose() * E0;
  const MATRIX3 H00 = term1 + term2;

  const std::tuple<MATRIX3, MATRIX3, MATRIX3> BH11 = binormalHessian11(E);
  term1.col(0) = std::get<0>(BH11).transpose() * rhs;
  term1.col(1) = std::get<1>(BH11).transpose() * rhs;
  term1.col(2) = std::get<2>(BH11).transpose() * rhs;
  //term1 *= 2.0;
  //term2 = 2.0 * E1.transpose() * M * B * M.transpose() * E1;
  term2 = E1.transpose() * M * B * M.transpose() * E1;
  const MATRIX3 H11 = term1 + term2;

  const std::tuple<MATRIX3, MATRIX3, MATRIX3> BH01 = binormalHessian01(E);
  term1.col(0) = std::get<0>(BH01).transpose() * rhs;
  term1.col(1) = std::get<1>(BH01).transpose() * rhs;
  term1.col(2) = std::get<2>(BH01).transpose() * rhs;
  //term1 *= 2.0;
  //term2 = 2.0 * E0.transpose() * M * B * M.transpose() * E1;
  term2 = E0.transpose() * M * B * M.transpose() * E1;
  const MATRIX3 H01 = term1 + term2;

  MATRIX6 H;
  H.block<3,3>(0,0) = H00;
  H.block<3,3>(3,0) = H01.transpose();
  H.block<3,3>(0,3) = H01;
  H.block<3,3>(3,3) = H11;

  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX2 SIN_BENDING::twistHessian(const MATRIX3x2& E, 
                                  const MATRIX2& B, 
                                  const MATRIX3x2& M0, 
                                  const MATRIX3x2& M1, 
                                  const VECTOR2& kappa, 
                                  const VECTOR2& kappaBar) const
{
  const VECTOR2 rhs = B * (kappa - kappaBar);
  const VECTOR3 ks = binormal(E);

  MATRIX2 T;
  T << 0, -1, 1, 0;

  MATRIX2 H;

  // Eigen really doesn't like doing this all in one shot
  H(0,0) = ks.transpose() * M0 * B * M0.transpose() * ks;
  H(0,0) *= 0.5;
  H(0,0) -= ((M0 * T).transpose() * ks).dot(rhs);

  H(1,1) = ks.transpose() * M1 * B * M1.transpose() * ks;
  H(1,1) *= 0.5;
  H(1,1) -= ((M1 * T).transpose() * ks).dot(rhs);

  H(0,1) = (ks.transpose() * M0 * B * M1.transpose() * ks);
  H(0,1) *= 0.5;

  H(1,0) = H(0,1);
  //return H;
  return 0.5 * H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<MATRIX3x2, MATRIX3x2> SIN_BENDING::mixedHessian(const MATRIX3x2& E, 
                                                          const MATRIX2& B, 
                                                          const MATRIX3x2& M0, 
                                                          const MATRIX3x2& M1, 
                                                          const MATRIX3x2& M, 
                                                          const VECTOR2& kappa, 
                                                          const VECTOR2& kappaBar) const
{
  const VECTOR2 rhs = B * (kappa - kappaBar);
  const std::pair<MATRIX3, MATRIX3> sinGradient = binormalGradient(E);
  const MATRIX3& E0 = sinGradient.first;
  const MATRIX3& E1 = sinGradient.second;

  const VECTOR3 ks = binormal(E);

  MATRIX3x2 H0, H1;
  H0.col(0) = -(E0.transpose() * M0 * rhs) - (E0.transpose() * M * B * M0.transpose() * ks);
  H0.col(1) = -(E0.transpose() * M1 * rhs) - (E0.transpose() * M * B * M1.transpose() * ks);

  H1.col(0) = -(E1.transpose() * M0 * rhs) - (E1.transpose() * M * B * M0.transpose() * ks);
  H1.col(1) = -(E1.transpose() * M1 * rhs) - (E1.transpose() * M * B * M1.transpose() * ks);

  H0 *= 0.5;
  H1 *= 0.5;
  return std::make_pair(H0,H1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR11 SIN_BENDING::gradient(const MATRIX3x2& E, 
                               const MATRIX2& B, 
                               const MATRIX3x2& M0, 
                               const MATRIX3x2& M1, 
                               const MATRIX3x2& M, 
                               const VECTOR2& kappa, 
                               const VECTOR2& kappaBar) const
{
  const MATRIX3x2 grad = edgeGradient(E, B, M, kappa, kappaBar);

  VECTOR9 perVertex = dedx().transpose() * flatten(grad);
  VECTOR2 perEdge = twistGradient(E, B, M0, M1, kappa, kappaBar);

  VECTOR11 result;
  result.setZero();
  result.segment<3>(0) = perVertex.segment<3>(0);
  result.segment<3>(4) = perVertex.segment<3>(3);
  result.segment<3>(8) = perVertex.segment<3>(6);

  result[3] = perEdge[0];
  result[7] = perEdge[1];

  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX11 SIN_BENDING::hessian(const MATRIX3x2& E, 
                              const MATRIX2& B, 
                              const MATRIX3x2& M0, 
                              const MATRIX3x2& M1, 
                              const MATRIX3x2& M, 
                              const VECTOR2& kappa, 
                              const VECTOR2& kappaBar) const
{
  MATRIX11 result;
  result.setZero();

  const MATRIX6 H = edgeHessian(E, B, M, kappaBar);
  const MATRIX6x9 changeOfBasis = dedx();

  // diagonal blocks
  result.block<9,9>(0,0) = changeOfBasis.transpose() * H * changeOfBasis;
  result.block<2,2>(9,9) = twistHessian(E, B, M0, M1, kappa, kappaBar);

  // off-diagonal terms
  std::pair<MATRIX3x2,MATRIX3x2> mixed = mixedHessian(E, B, M0, M1, M, kappa, kappaBar);
  const MATRIX3x2& H0 = mixed.first;
  const MATRIX3x2& H1 = mixed.second;

  VECTOR6 E0;
  E0.segment<3>(0) = H0.col(0);
  E0.segment<3>(3) = H1.col(0);
  
  VECTOR6 E1;
  E1.segment<3>(0) = H0.col(1);
  E1.segment<3>(3) = H1.col(1);

  const VECTOR9 X0 = changeOfBasis.transpose() * E0;
  const VECTOR9 X1 = changeOfBasis.transpose() * E1;

  result.block<9,1>(0,9)  = X0;
  result.block<9,1>(0,10) = X1;
  
  result.block<1,9>(9,0)  = X0.transpose();
  result.block<1,9>(10,0) = X1.transpose();

  // need to rewrite this if the energy starts looking promising
  MATRIX11 thetasLast = result;
  result.col(0)  = thetasLast.col(0);
  result.col(1)  = thetasLast.col(1);
  result.col(2)  = thetasLast.col(2);
  result.col(3)  = thetasLast.col(9);
  result.col(4)  = thetasLast.col(3);
  result.col(5)  = thetasLast.col(4);
  result.col(6)  = thetasLast.col(5);
  result.col(7)  = thetasLast.col(10);
  result.col(8)  = thetasLast.col(6);
  result.col(9)  = thetasLast.col(7);
  result.col(10) = thetasLast.col(8);
  
  thetasLast = result;
  result.row(0)  = thetasLast.row(0);
  result.row(1)  = thetasLast.row(1);
  result.row(2)  = thetasLast.row(2);
  result.row(3)  = thetasLast.row(9);
  result.row(4)  = thetasLast.row(3);
  result.row(5)  = thetasLast.row(4);
  result.row(6)  = thetasLast.row(5);
  result.row(7)  = thetasLast.row(10);
  result.row(8)  = thetasLast.row(6);
  result.row(9)  = thetasLast.row(7);
  result.row(10) = thetasLast.row(8);

  return result;
}

} // STRAND 
} // ANGLE
