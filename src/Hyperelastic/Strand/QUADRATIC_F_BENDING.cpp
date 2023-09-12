#include "QUADRATIC_F_BENDING.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

#define USING_GAUSS_NEWTON_BENDING_STRAND 0

QUADRATIC_F_BENDING::QUADRATIC_F_BENDING(const REAL& mu) :
  TAN_BENDING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string QUADRATIC_F_BENDING::name() const
{ 
  return std::string("QUADRATIC_F_BENDING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_F_BENDING::psi(const MATRIX3x2& E, const REAL& theta0) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL theta = acos(cosTheta);
  const REAL diff = theta - theta0;

  return 0.5 * _mu * diff * diff;
  //return 0.5 * _mu * theta;
}


REAL QUADRATIC_F_BENDING::psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  // this mirrors computation in ISOTROPIC_THETA
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  REAL theta = acos(cosTheta);
  theta = inverted? -theta:theta;
  REAL sign = inverted?-1.0:1.0;
  const REAL diff = theta - theta0;

  return 0.5 * _mu * diff * diff;
  //return 0.5 * _mu * theta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// gradient of theta
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_F_BENDING::gradient(const MATRIX3x2& E) const
{
	const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const VECTOR3 crs = e0.cross(e1);
  const VECTOR3 n = crs.normalized();
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  const VECTOR3 g0 = e0perp / e0.dot(e0);
  const VECTOR3 g1 = -e1perp / e1.dot(e1);
  MATRIX3x2 G; G << g0, g1;
  return G;
}

static REAL clamp(const REAL& input, const REAL& bottom, const REAL& top)
{
  REAL result = (input > bottom) ? input : bottom;
  return result < top ? result : top;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// returns Eijk * v
///////////////////////////////////////////////////////////////////////////////////////////////////
// MATRIX3 getCrossMatrix(const VECTOR3 v)
// {
//   MATRIX3 M;
//   M << 0, v[2], -v[1],
//        -v[2], 0, v[0],
//        v[1], -v[0], 0;
//   return M;
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 QUADRATIC_F_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  //const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);
  // cout<<"theta: "<<theta<<" theta0: "<<theta0<<endl;
  const MATRIX3x2 Pfinal = _mu * (theta - theta0) * gradient(E);

  return Pfinal;
}

MATRIX3x2 QUADRATIC_F_BENDING::PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  //const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  REAL theta = acos(cosTheta);
  theta = inverted? -theta:theta;
  REAL sign = inverted?-1.0:1.0;
  const MATRIX3x2 Pfinal = sign * _mu * (theta - theta0) * gradient(E);

  return Pfinal;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_F_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0) const
{
	// compute the Hessian of theta
	const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const VECTOR3 crs = e0.cross(e1);
  const REAL cnorm = crs.norm();
  const VECTOR3 n = crs.normalized();
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  const MATRIX3 N = getCrossMatrix(n);
  const MATRIX3 E0 = getCrossMatrix(e0);
  const MATRIX3 E1 = getCrossMatrix(e1);
  const REAL dot0 = e0.dot(e0);
  const REAL dot1 = e1.dot(e1);

  // d Psi d theta
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);
  const REAL g = _mu * (theta - theta0);
  const REAL h = _mu;

  // the gradient of theta
  VECTOR6 G;
  const VECTOR3 g0 = e0perp / dot0, g1 = -e1perp / dot1;
  G << g0, g1;

  #if USING_GAUSS_NEWTON_BENDING_STRAND
    // cout<<"--------gauss newton.\n";
    return h * G * G.transpose();
  #endif

  // derivatives of n
  MATRIX3 pnpe0, pnpe1;
  if(cnorm < 1e-16){
    pnpe0.setZero();
    pnpe1.setZero();
  }
  else{
    pnpe0 = (E1 * cnorm - crs * e1perp.transpose())/cnorm/cnorm;
    pnpe1 = (- E0 * cnorm + crs * e0perp.transpose())/cnorm/cnorm;
  }
  // compute the blocks
  const MATRIX3 B00 = (dot0 * N - 2 * e0perp * e0.transpose() + (-E0 * pnpe0) * dot0)/dot0/dot0;
  const MATRIX3 B01 = -E0 * pnpe1 / dot0;
  const MATRIX3 B11 = (-dot1 * N + 2 * e1perp * e1.transpose() + (E1 * pnpe1) * dot1)/dot1/dot1;
  const MATRIX3 B10 = E1 * pnpe0 / dot1;
  MATRIX6 hTehta;
  hTehta << B00, B01,
  					B10, B11;



  return h * G * G.transpose() + g * hTehta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_F_BENDING::hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  // compute the Hessian of theta
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const VECTOR3 crs = e0.cross(e1);
  const REAL cnorm = crs.norm();
  const VECTOR3 n = crs.normalized();
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  const MATRIX3 N = getCrossMatrix(n);
  const MATRIX3 E0 = getCrossMatrix(e0);
  const MATRIX3 E1 = getCrossMatrix(e1);
  const REAL dot0 = e0.dot(e0);
  const REAL dot1 = e1.dot(e1);
  // derivatives of n
  MATRIX3 pnpe0, pnpe1;
  if(cnorm < 1e-16){
    pnpe0.setZero();
    pnpe1.setZero();
  }
  else{
    pnpe0 = (E1 * cnorm - crs * e1perp.transpose())/cnorm/cnorm;
    pnpe1 = (- E0 * cnorm + crs * e0perp.transpose())/cnorm/cnorm;
  }
  // compute the blocks
  const MATRIX3 B00 = (dot0 * N - 2 * e0perp * e0.transpose() + (-E0 * pnpe0) * dot0)/dot0/dot0;
  const MATRIX3 B01 = -E0 * pnpe1 / dot0;
  const MATRIX3 B11 = (-dot1 * N + 2 * e1perp * e1.transpose() + (E1 * pnpe1) * dot1)/dot1/dot1;
  const MATRIX3 B10 = E1 * pnpe0 / dot1;
  MATRIX6 hTehta;
  hTehta << B00, B01,
            B10, B11;

  // the gradient of theta
  VECTOR6 G;
  const VECTOR3 g0 = e0perp / dot0, g1 = -e1perp / dot1;
  G << g0, g1;

  // d Psi d theta
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  REAL theta = acos(cosTheta);
  theta = inverted? -theta:theta;
  REAL sign = inverted?-1.0:1.0;
  const REAL g = _mu * (theta - theta0);
  const REAL h = _mu;
  #if USING_GAUSS_NEWTON_BENDING_STRAND
    // cout<<"--------gauss newton.\n";
    return h * G * G.transpose();
  #endif
  return h * G * G.transpose() + sign * g * hTehta;
}

MATRIX6 QUADRATIC_F_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  const REAL theta = acos(cosTheta);
  const VECTOR3 crs = e0.cross(e1);
  REAL cnorm = crs.norm();
  VECTOR3 n = crs.normalized();
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  REAL dot0 = e0.dot(e0), dot1 = e1.dot(e1);
  REAL norm0 = e0.norm(), norm1 = e1.norm();

  const REAL cosine = cos(theta * 0.5), sine = sin(theta * 0.5), tangent = tan(theta * 0.5);
  // partial psi partial theta
  const REAL g = _mu * (theta - theta0);
  const REAL h = _mu;

  // cout << "g: "<<g<<"h: "<<h<<endl;

  // the two eigenvectors composed of n
  const REAL beta = (norm1/norm0 - norm0/norm1) * cos(theta);
  const REAL alpha = (- beta + sqrt(beta * beta + 4.0)) * 0.5;
  VECTOR6 q4, q5;
  REAL lambda4, lambda5;
  q4 << alpha * n, n;
  q5 << n, -alpha * n;
  if(abs(theta) < 1e-16){
    lambda4 = 0.0;
    lambda5 = 0.0;
  }
  else{
    lambda4 = g*(1.0/dot1/tan(theta) - alpha/norm0/norm1/sin(theta));
    lambda5 = g*(1.0/dot1/tan(theta) + 1.0/alpha/norm0/norm1/sin(theta));
  }
  // the rest 4
  const REAL n0 = h / dot0, n1 = h / dot1, w = g / h;
  const REAL w2 = w * w, nDiff = n0 - n1, nSum = n0 + n1;
  REAL gamma = nDiff/nSum; gamma = gamma * gamma;
  const REAL r = sqrt(4.0 * w2 * gamma + 1.0) * nSum;
  const REAL Rminus = sqrt(2.0 * (2.0 * w2 + 1.0 - r/nSum)) * nSum;
  const REAL Rplus = sqrt(2.0 * (2.0 * w2 + 1.0 + r/nSum)) * nSum;
// cout << "r: "<<r<<"Rminus: "<<Rminus<<"nSum: "<<nSum<<endl;
  const REAL lambda0 = 0.25 * (-r - Rminus + nSum);
  const REAL lambda1 = 0.25 * (-r + Rminus + nSum);
  const REAL lambda2 = 0.25 * (r - Rplus + nSum);
  const REAL lambda3 = 0.25 * (r + Rplus + nSum);

  VECTOR4 q0, q1, q2, q3;
  // q0.setZero(); q1.setZero(); q2.setZero(); q3.setZero(); 
  q0[0] = w/n1*(nDiff * nDiff - nSum * r - nDiff * Rminus);
  q0[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r+0.5*(nDiff-r)*Rminus)/n1;
  q0[2] = 4.0 * n1 * w;
  q0[3] = nSum - r - Rminus;
  q1[0] = w/n1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/n1;
  q1[2] = 4.0 * n1 * w;
  q1[3] = nSum - r + Rminus;
  q2[0] = w/n1*(nDiff * nDiff + nSum * r - nDiff * Rplus);
  q2[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r+0.5*(nDiff+r)*Rplus)/n1;
  q2[2] = 4.0 * n1 * w;
  q2[3] = nSum + r - Rplus;
  q3[0] = w/n1*(nDiff * nDiff + nSum * r + nDiff * Rplus);
  q3[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r-0.5*(nDiff+r)*Rplus)/n1;
  q3[2] = 4.0 * n1 * w;
  q3[3] = nSum + r + Rplus;

  MATRIX6x4 A;
  A << e0, e0perp, VECTOR3::Zero(), VECTOR3::Zero(),
        VECTOR3::Zero(), VECTOR3::Zero(), e1, e1perp;
  VECTOR6 qq0 = A * q0, qq1 = A * q1, qq2 = A * q2, qq3 = A * q3;
  qq0 = qq0.normalized(); qq1 = qq1.normalized(); qq2 = qq2.normalized(); 
  qq3 = qq3.normalized(); q4 = q4.normalized(); q5 = q5.normalized();
  MATRIX6 Q; Q << qq0, qq1, qq2, qq3, q4, q5;
  VECTOR6 lambda; lambda << lambda0, lambda1, lambda2, lambda3, lambda4, lambda5;
  // cout<<"lambda: "<<lambda << endl;
  for (unsigned int x = 0; x < 6; x++)
    lambda[x] = (lambda[x] > 0.0) ? lambda[x] : 0.0;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();
  return H;
}

MATRIX6 QUADRATIC_F_BENDING::clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const
{
  const VECTOR3 e0 = E.col(0);
  const VECTOR3 e1 = E.col(1);
  const REAL cosTheta = clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0);
  REAL theta = acos(cosTheta);
  theta = inverted? -theta:theta;
  REAL sign = inverted?-1.0:1.0;
  const VECTOR3 crs = e0.cross(e1);
  REAL cnorm = crs.norm();
  VECTOR3 n = crs.normalized();
  const VECTOR3 e0perp = e0.cross(n);
  const VECTOR3 e1perp = e1.cross(n);
  REAL dot0 = e0.dot(e0), dot1 = e1.dot(e1);
  REAL norm0 = e0.norm(), norm1 = e1.norm();

  const REAL cosine = cos(theta * 0.5), sine = sin(theta * 0.5), tangent = tan(theta * 0.5);
  // partial psi partial theta
  const REAL g = _mu * (theta - theta0);
  const REAL h = _mu;

  // cout << "g: "<<g<<"h: "<<h<<endl;

  // the two eigenvectors composed of n
  const REAL beta = (norm1/norm0 - norm0/norm1) * cos(theta);
  const REAL alpha = (- beta + sqrt(beta * beta + 4.0)) * 0.5;
  VECTOR6 q4, q5;
  REAL lambda4, lambda5;
  q4 << alpha * n, n;
  q5 << n, -alpha * n;
  if(abs(theta) < 1e-16){
    lambda4 = 0.0;
    lambda5 = 0.0;
  }
  else{
    lambda4 = sign * g*(1.0/dot1/tan(theta) - alpha/norm0/norm1/sin(theta));
    lambda5 = sign * g*(1.0/dot1/tan(theta) + 1.0/alpha/norm0/norm1/sin(theta));
  }
  // the rest 4
  const REAL n0 = h / dot0, n1 = h / dot1, w = sign * g / h;
  const REAL w2 = w * w, nDiff = n0 - n1, nSum = n0 + n1;
  REAL gamma = nDiff/nSum; gamma = gamma * gamma;
  const REAL r = sqrt(4.0 * w2 * gamma + 1.0) * nSum;
  const REAL Rminus = sqrt(2.0 * (2.0 * w2 + 1.0 - r/nSum)) * nSum;
  const REAL Rplus = sqrt(2.0 * (2.0 * w2 + 1.0 + r/nSum)) * nSum;
// cout << "r: "<<r<<"Rminus: "<<Rminus<<"nSum: "<<nSum<<endl;
  const REAL lambda0 = 0.25 * (-r - Rminus + nSum);
  const REAL lambda1 = 0.25 * (-r + Rminus + nSum);
  const REAL lambda2 = 0.25 * (r - Rplus + nSum);
  const REAL lambda3 = 0.25 * (r + Rplus + nSum);

  VECTOR4 q0, q1, q2, q3;
  // q0.setZero(); q1.setZero(); q2.setZero(); q3.setZero(); 
  q0[0] = w/n1*(nDiff * nDiff - nSum * r - nDiff * Rminus);
  q0[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r+0.5*(nDiff-r)*Rminus)/n1;
  q0[2] = 4.0 * n1 * w;
  q0[3] = nSum - r - Rminus;
  q1[0] = w/n1*(nDiff * nDiff - nSum * r + nDiff * Rminus);
  q1[1] = (nSum*(-n0-2.0*w2*nDiff)+n0*r-0.5*(nDiff-r)*Rminus)/n1;
  q1[2] = 4.0 * n1 * w;
  q1[3] = nSum - r + Rminus;
  q2[0] = w/n1*(nDiff * nDiff + nSum * r - nDiff * Rplus);
  q2[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r+0.5*(nDiff+r)*Rplus)/n1;
  q2[2] = 4.0 * n1 * w;
  q2[3] = nSum + r - Rplus;
  q3[0] = w/n1*(nDiff * nDiff + nSum * r + nDiff * Rplus);
  q3[1] = (nSum*(-n0-2.0*w2*nDiff)-n0*r-0.5*(nDiff+r)*Rplus)/n1;
  q3[2] = 4.0 * n1 * w;
  q3[3] = nSum + r + Rplus;

  MATRIX6x4 A;
  A << e0, e0perp, VECTOR3::Zero(), VECTOR3::Zero(),
        VECTOR3::Zero(), VECTOR3::Zero(), e1, e1perp;
  VECTOR6 qq0 = A * q0, qq1 = A * q1, qq2 = A * q2, qq3 = A * q3;
  qq0 = qq0.normalized(); qq1 = qq1.normalized(); qq2 = qq2.normalized(); 
  qq3 = qq3.normalized(); q4 = q4.normalized(); q5 = q5.normalized();
  MATRIX6 Q; Q << qq0, qq1, qq2, qq3, q4, q5;
  VECTOR6 lambda; lambda << lambda0, lambda1, lambda2, lambda3, lambda4, lambda5;
  // cout<<"lambda: "<<lambda << endl;
  for (unsigned int x = 0; x < 6; x++)
    lambda[x] = (lambda[x] > 0.0) ? lambda[x] : 0.0;

  const MATRIX6 H = Q * lambda.asDiagonal() * Q.transpose();
  return H;
}



} // STRAND 
} // ANGLE