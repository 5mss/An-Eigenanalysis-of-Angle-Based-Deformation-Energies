using namespace HOBAK;
using namespace std;

const REAL theta0 = M_PI * 0.5;

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const HOBAK::STRAND::STRETCHING* material, const VECTOR3& f)
{
  using namespace HOBAK;
  using namespace std;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3 dPdF = material->hessian(f);
  VECTOR3 p = material->PK1(f);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX3 finiteDiff;

    // for each of the degrees of the freedom
    for (int x = 0; x < 3; x++)
    {
      VECTOR3 fnew = f;
      fnew[x] += eps;

      // get the new psi
      VECTOR3 pnew = material->PK1(fnew);

      // store the finite difference
      VECTOR3 diff = (pnew - p) / eps;
      finiteDiff.col(x) = diff;
    }

    MATRIX3 diff = dPdF - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / p.norm())) / 9.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    MATRIX3 div = finiteDiff;
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
        div(x,y) = div(x,y) / dPdF(x,y);
    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " dPdF: " << endl << dPdF << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      cout << " div: " << endl << div << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestSpatialHessian(const HOBAK::STRAND::STRETCHING* material, const vector<VECTOR3>& p)
{
  using namespace HOBAK;
  using namespace std;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand spatial Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  const REAL dmInv = 1.2345;
  const MATRIX6 H = material->spatialHessian(p, dmInv);
  const VECTOR6 g = material->spatialGradient(p, dmInv);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX6 finiteDiff;

    // for each of the degrees of the freedom
    int index = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++, index++)
      {
        vector<VECTOR3> pNew = p;
        pNew[y][x] += eps;

        // get the new psi
        VECTOR6 gNew = material->spatialGradient(pNew, dmInv);

        // store the finite difference
        VECTOR6 diff = (gNew - g) / eps;
        finiteDiff.col(index) = diff;
      }

    MATRIX6 diff = H - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / H.norm())) / 36.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    MATRIX6 div = finiteDiff;
    for (int y = 0; y < 6; y++)
      for (int x = 0; x < 6; x++)
        div(x,y) = div(x,y) / H(x,y);
    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      cout << " div: " << endl << div << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////////////////////////
// Third vector to complete a orthonormal frame
///////////////////////////////////////////////////////////////////////
static VECTOR3 curvatureBinormal(const VECTOR3& e0,
                                 const VECTOR3& e1)
{
  const REAL denom = e0.norm() *  e1.norm() + e0.dot(e1);
  return (2.0 / denom) * e0.cross(e1);
}

#if 0
///////////////////////////////////////////////////////////////////////
// Panetta's approach to parallel transport
///////////////////////////////////////////////////////////////////////
static VECTOR3 parallelTransport(const VECTOR3& t0,
                                 const VECTOR3& t1,
                                 const VECTOR3& v)
{
  const VECTOR3 sinTheta = t0.cross(t1);
  const REAL cosTheta = t0.dot(t1);
  const REAL onePlus = 1.0 + cosTheta;

  // if it's the same vector, but pointing in the opposite direction,
  // the frame didn't change (didn't rotate 180?)
  if (fabs(onePlus) < 1e-14)
    return v;

  // if it's the same vector, then the frame didn't change
  if ((t0 - t1).cwiseAbs().maxCoeff() == 0) return v;

  return (sinTheta.dot(v) / onePlus) * sinTheta + sinTheta.cross(v) + cosTheta * v;
}
#endif

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the spatial gradient of sin-based bending
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestThetaGradientSinBending(const HOBAK::STRAND::SIN_BENDING& bending)
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand twist gradient for " << bending.name().c_str() << " bending" << endl;
  cout << "=============================================================== " << endl;

  VECTOR9 positions;
  VECTOR3 xBar0(0,0,0);
  VECTOR3 xBar1(1,1,0);
  VECTOR3 xBar2(2,0,1);

  positions.block<3,1>(0,0) = xBar0;
  positions.block<3,1>(3,0) = xBar1;
  positions.block<3,1>(6,0) = xBar2;

  VECTOR2 thetas;
  thetas[0] =  M_PI * 0.25;
  thetas[1] = -M_PI * 0.25;

  VECTOR3 edgeBar0 = xBar1 - xBar0;
  VECTOR3 edgeBar1 = xBar2 - xBar1;
 
  // TODO: noise this up
  const VECTOR3 x0 = randomVector3();  
  const VECTOR3 x1 = randomVector3();  
  const VECTOR3 x2 = randomVector3();  
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;

  MATRIX3x2 Ebar;
  Ebar.col(0) = edgeBar0;
  Ebar.col(1) = edgeBar1;

  MATRIX3x2 E;
  E.col(0) = edge0;
  E.col(1) = edge1;
 
  const VECTOR3 kbBar = bending.binormal(edgeBar0, edgeBar1);
  const VECTOR3 kb    = bending.binormal(edge0, edge1);

  const VECTOR3 d1_0(0,1,0);
  const VECTOR3 d2_0(0,0,1);
  const VECTOR3 d1_1(0,1,0);
  const VECTOR3 d2_1(0,0,1);

  const VECTOR3 m1_0 =  cos(thetas[0]) * d1_0 + sin(thetas[0]) * d2_0;
  const VECTOR3 m2_0 = -sin(thetas[0]) * d1_0 + cos(thetas[0]) * d2_0;
  const VECTOR3 m1_1 =  cos(thetas[1]) * d1_0 + sin(thetas[1]) * d2_0;
  const VECTOR3 m2_1 = -sin(thetas[1]) * d1_0 + cos(thetas[1]) * d2_0;
  
  MATRIX3x2 M;
  M.col(0) =  0.5 * (m2_0 + m2_1);
  M.col(1) = -0.5 * (m1_0 + m1_1);

  const VECTOR2 kappaBar = M.transpose() * kbBar;
  const VECTOR2 kappa = M.transpose() * kb;

  MATRIX2 B;
  B.setZero();
  B(0,0) = 2.0;
  B(1,1) = 1.0;

  MATRIX3x2 M0;
  M0.col(0) = m1_0;
  M0.col(1) = m2_0;
    
  MATRIX3x2 M1;
  M1.col(0) = m1_1;
  M1.col(1) = m2_1;

  const REAL psi0 = bending.psi(B, kappa, kappaBar);
  const VECTOR2 g = bending.twistGradient(E, B, M0, M1, kappa, kappaBar);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  VECTOR2 finiteDiff;
  while (eps > 1e-8)
  {
    finiteDiff.setZero();

    // for each of the degrees of the freedom
    for (int i = 0; i < 2; i++)
    {
      // add the pertubation
      VECTOR2 thetasNew = thetas;
      thetasNew[i] += eps;

      const VECTOR3 m1_0new =  cos(thetasNew[0]) * d1_0 + sin(thetasNew[0]) * d2_0;
      const VECTOR3 m2_0new = -sin(thetasNew[0]) * d1_0 + cos(thetasNew[0]) * d2_0;
      const VECTOR3 m1_1new =  cos(thetasNew[1]) * d1_0 + sin(thetasNew[1]) * d2_0;
      const VECTOR3 m2_1new = -sin(thetasNew[1]) * d1_0 + cos(thetasNew[1]) * d2_0;
      MATRIX3x2 Mnew;
      Mnew.col(0) =  0.5 * (m2_0new + m2_1new);
      Mnew.col(1) = -0.5 * (m1_0new + m1_1new);
      const VECTOR2 kappaNew = Mnew.transpose() * kb;
      
      // get the new psi
      double psi = bending.psi(B, kappaNew, kappaBar);

      // store the finite difference
      finiteDiff[i] = (psi - psi0) / eps;
    }

    VECTOR2 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 2.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 2.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " g: " << endl << g << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
    cout << " analytic: " << endl << g << endl;
    cout << " numerical: " << endl << finiteDiff << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the edge gradient of sin-based bending
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeGradientSinBending(const HOBAK::STRAND::SIN_BENDING& bending)
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand edge gradient for " << bending.name().c_str() << " bending" << endl;
  cout << "=============================================================== " << endl;

  VECTOR9 positions;
  VECTOR3 xBar0(0,0,0);
  VECTOR3 xBar1(1,1,0);
  VECTOR3 xBar2(2,0,1);

  positions.block<3,1>(0,0) = xBar0;
  positions.block<3,1>(3,0) = xBar1;
  positions.block<3,1>(6,0) = xBar2;

  VECTOR2 thetas;
  thetas[0] =  M_PI * 0.25;
  thetas[1] = -M_PI * 0.25;
  //thetas[0] = 0;
  //thetas[1] = 0;

  VECTOR3 edgeBar0 = xBar1 - xBar0;
  VECTOR3 edgeBar1 = xBar2 - xBar1;

  const VECTOR3 x0 = randomVector3();  
  const VECTOR3 x1 = randomVector3();  
  const VECTOR3 x2 = randomVector3();  
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;

  MATRIX3x2 Ebar;
  Ebar.col(0) = edgeBar0;
  Ebar.col(1) = edgeBar1;

  MATRIX3x2 E;
  E.col(0) = edge0;
  E.col(1) = edge1;
 
  const VECTOR3 kbBar = bending.binormal(edgeBar0, edgeBar1);
  const VECTOR3 kb    = bending.binormal(edge0, edge1);

  const VECTOR3 d1_0(0,1,0);
  const VECTOR3 d2_0(0,0,1);
  const VECTOR3 d1_1(0,1,0);
  const VECTOR3 d2_1(0,0,1);

  const VECTOR3 m1_0 =  cos(thetas[0]) * d1_0 + sin(thetas[0]) * d2_0;
  const VECTOR3 m2_0 = -sin(thetas[0]) * d1_0 + cos(thetas[0]) * d2_0;
  const VECTOR3 m1_1 =  cos(thetas[1]) * d1_0 + sin(thetas[1]) * d2_0;
  const VECTOR3 m2_1 = -sin(thetas[1]) * d1_0 + cos(thetas[1]) * d2_0;
  
  MATRIX3x2 M;
  M.col(0) =  0.5 * (m2_0 + m2_1);
  M.col(1) = -0.5 * (m1_0 + m1_1);

  const VECTOR2 kappaBar = M.transpose() * kbBar;
  const VECTOR2 kappa = M.transpose() * kb;

  MATRIX2 B;
  B.setZero();
  B(0,0) = 2.0;
  B(1,1) = 1.0;

  const REAL psi0 = bending.psi(B, kappa, kappaBar);
  const MATRIX3x2 g = bending.edgeGradient(E, B, M, kappa, kappaBar);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  MATRIX3x2 finiteDiff;
  while (eps > 1e-8)
  {
    finiteDiff.setZero();

    // for each of the degrees of the freedom
    int index = 0;
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 3; i++, index++)
      {
        // add the pertubation
        MATRIX3x2 Enew = E;
        Enew(i,j) += eps;
  
        const VECTOR3 kbNew = bending.binormal(Enew);
        const VECTOR2 kappaNew = M.transpose() * kbNew;
        
        // get the new psi
        double psi = bending.psi(B, kappaNew, kappaBar);

        // store the finite difference
        finiteDiff(i,j) = (psi - psi0) / eps;
      }

    MATRIX3x2 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 6.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 6.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " g: " << endl << g << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
    //cout << " analytic: " << endl << g << endl;
    //cout << " numerical: " << endl << finiteDiff << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the edge hessian of sin-based bending
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessianSinBending(const HOBAK::STRAND::SIN_BENDING& bending)
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand hessian for " << bending.name().c_str() << " bending" << endl;
  cout << "=============================================================== " << endl;

  VECTOR3 xBar0(0,0,0);
  VECTOR3 xBar1(1,1,0);
  VECTOR3 xBar2(2,0,1);

  VECTOR9 positionsBar;
  positionsBar.block<3,1>(0,0) = xBar0;
  positionsBar.block<3,1>(3,0) = xBar1;
  positionsBar.block<3,1>(6,0) = xBar2;

  VECTOR2 thetas;
  thetas[0] =  M_PI * 0.25;
  thetas[1] = -M_PI * 0.25;
  //thetas[0] = 0;
  //thetas[1] = 0;

  VECTOR3 edgeBar0 = xBar1 - xBar0;
  VECTOR3 edgeBar1 = xBar2 - xBar1;

  const VECTOR3 x0 = randomVector3();  
  const VECTOR3 x1 = randomVector3();  
  const VECTOR3 x2 = randomVector3();  
  //const VECTOR3 x0 = VECTOR3(0,0,0);
  //const VECTOR3 x1 = VECTOR3(1,0,0);
  //const VECTOR3 x2 = VECTOR3(2,0,0);
  VECTOR9 positions;
  positions.block<3,1>(0,0) = x0;
  positions.block<3,1>(3,0) = x1;
  positions.block<3,1>(6,0) = x2;
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;

  MATRIX3x2 Ebar;
  Ebar.col(0) = edgeBar0;
  Ebar.col(1) = edgeBar1;

  MATRIX3x2 E;
  E.col(0) = edge0;
  E.col(1) = edge1;
  
  const VECTOR3 kbBar = bending.binormal(edgeBar0, edgeBar1);
  const VECTOR3 kb    = bending.binormal(edge0, edge1);

  const VECTOR3 d1_0(0,1,0);
  const VECTOR3 d2_0(0,0,1);
  const VECTOR3 d1_1(0,1,0);
  const VECTOR3 d2_1(0,0,1);

  const VECTOR3 m1_0 =  cos(thetas[0]) * d1_0 + sin(thetas[0]) * d2_0;
  const VECTOR3 m2_0 = -sin(thetas[0]) * d1_0 + cos(thetas[0]) * d2_0;
  const VECTOR3 m1_1 =  cos(thetas[1]) * d1_0 + sin(thetas[1]) * d2_0;
  const VECTOR3 m2_1 = -sin(thetas[1]) * d1_0 + cos(thetas[1]) * d2_0;
  
  MATRIX3x2 M;
  M.col(0) =  0.5 * (m2_0 + m2_1);
  M.col(1) = -0.5 * (m1_0 + m1_1);

  const VECTOR2 kappaBar = M.transpose() * kbBar;
  const VECTOR2 kappa = M.transpose() * kb;

  MATRIX2 B;
  B.setZero();
  B(0,0) = 2.0;
  B(1,1) = 1.0;

  MATRIX3x2 M0;
  M0.col(0) = m1_0;
  M0.col(1) = m2_0;
    
  MATRIX3x2 M1;
  M1.col(0) = m1_1;
  M1.col(1) = m2_1;

  const VECTOR11 g0 = bending.gradient(E, B, M0, M1, M, kappa, kappaBar);
  const MATRIX11 H = bending.hessian(E, B, M0, M1, M, kappa, kappaBar);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  MATRIX11 finiteDiff;
  while (eps > 1e-8)
  {
    finiteDiff.setZero();

    // position degrees
    for (int i = 0; i < 11; i++)
    {
      if (i != 3 && i != 7)
      {
        int index = i;
        if (i > 3) index--;
        if (i > 7) index--;

        VECTOR9 positionsNew = positions;
        positionsNew[index] += eps;
    
        const VECTOR3 x0New = positionsNew.segment<3>(0);
        const VECTOR3 x1New = positionsNew.segment<3>(3);
        const VECTOR3 x2New = positionsNew.segment<3>(6);

        const VECTOR3 edge0New = x1New - x0New;
        const VECTOR3 edge1New = x2New - x1New;
        MATRIX3x2 Enew;
        Enew.col(0) = edge0New;
        Enew.col(1) = edge1New;

        const VECTOR3 kbNew = bending.binormal(Enew);
        const VECTOR2 kappaNew = M.transpose() * kbNew;
        
        // get the new psi
        const VECTOR11 g = bending.gradient(Enew, B, M0, M1, M, kappaNew, kappaBar);

        // store the finite difference
        finiteDiff.col(i) = (g - g0) / eps;
      }
      else
      {
        int index = 0;
        if (i == 7) index = 1;
        // add the pertubation
        VECTOR2 thetasNew = thetas;
        thetasNew[index] += eps;

        const VECTOR3 m1_0new =  cos(thetasNew[0]) * d1_0 + sin(thetasNew[0]) * d2_0;
        const VECTOR3 m2_0new = -sin(thetasNew[0]) * d1_0 + cos(thetasNew[0]) * d2_0;
        const VECTOR3 m1_1new =  cos(thetasNew[1]) * d1_0 + sin(thetasNew[1]) * d2_0;
        const VECTOR3 m2_1new = -sin(thetasNew[1]) * d1_0 + cos(thetasNew[1]) * d2_0;
        MATRIX3x2 Mnew;
        Mnew.col(0) =  0.5 * (m2_0new + m2_1new);
        Mnew.col(1) = -0.5 * (m1_0new + m1_1new);

        MATRIX3x2 M0new;
        M0new.col(0) = m1_0new;
        M0new.col(1) = m2_0new;
      
        MATRIX3x2 M1new;
        M1new.col(0) = m1_1new;
        M1new.col(1) = m2_1new;

        const VECTOR2 kappaNew = Mnew.transpose() * kb;
        
        // get the new psi
        const VECTOR11 g = bending.gradient(E, B, M0new, M1new, Mnew, kappaNew, kappaBar);

        // store the finite difference
        finiteDiff.col(i) = (g - g0) / eps;
      }
    }

    /*    
    // position degrees
    for (int i = 0; i < 9; i++)
    {
      VECTOR9 positionsNew = positions;
      positionsNew[i] += eps;
  
      const VECTOR3 x0New = positionsNew.segment<3>(0);
      const VECTOR3 x1New = positionsNew.segment<3>(3);
      const VECTOR3 x2New = positionsNew.segment<3>(6);

      const VECTOR3 edge0New = x1New - x0New;
      const VECTOR3 edge1New = x2New - x1New;
      MATRIX3x2 Enew;
      Enew.col(0) = edge0New;
      Enew.col(1) = edge1New;

      const VECTOR3 kbNew = bending.binormal(Enew);
      const VECTOR2 kappaNew = M.transpose() * kbNew;
      
      // get the new psi
      const VECTOR11 g = bending.gradient(Enew, B, M0, M1, M, kappaNew, kappaBar);

      // store the finite difference
      finiteDiff.col(i) = (g - g0) / eps;
    }
    for (int i = 0; i < 2; i++)
    {
      // add the pertubation
      VECTOR2 thetasNew = thetas;
      thetasNew[i] += eps;

      const VECTOR3 m1_0new =  cos(thetasNew[0]) * d1_0 + sin(thetasNew[0]) * d2_0;
      const VECTOR3 m2_0new = -sin(thetasNew[0]) * d1_0 + cos(thetasNew[0]) * d2_0;
      const VECTOR3 m1_1new =  cos(thetasNew[1]) * d1_0 + sin(thetasNew[1]) * d2_0;
      const VECTOR3 m2_1new = -sin(thetasNew[1]) * d1_0 + cos(thetasNew[1]) * d2_0;
      MATRIX3x2 Mnew;
      Mnew.col(0) =  0.5 * (m2_0new + m2_1new);
      Mnew.col(1) = -0.5 * (m1_0new + m1_1new);

      MATRIX3x2 M0new;
      M0new.col(0) = m1_0new;
      M0new.col(1) = m2_0new;
    
      MATRIX3x2 M1new;
      M1new.col(0) = m1_1new;
      M1new.col(1) = m2_1new;

      const VECTOR2 kappaNew = Mnew.transpose() * kb;
      
      // get the new psi
      const VECTOR11 g = bending.gradient(E, B, M0new, M1new, Mnew, kappaNew, kappaBar);

      // store the finite difference
      finiteDiff.col(9 + i) = (g - g0) / eps;
    }
    */

    MATRIX11 diff = H - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / H.norm())) / 121.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 121.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
    //cout << " analytic: " << endl << g << endl;
    //cout << " numerical: " << endl << finiteDiff << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the edge hessian of sin-based bending
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestEdgeHessianSinBending(const HOBAK::STRAND::SIN_BENDING& bending)
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand edge hessian for " << bending.name().c_str() << " bending" << endl;
  cout << "=============================================================== " << endl;

  VECTOR9 positions;
  VECTOR3 xBar0(0,0,0);
  VECTOR3 xBar1(1,1,0);
  VECTOR3 xBar2(2,0,1);

  positions.block<3,1>(0,0) = xBar0;
  positions.block<3,1>(3,0) = xBar1;
  positions.block<3,1>(6,0) = xBar2;

  VECTOR2 thetas;
  thetas[0] =  M_PI * 0.25;
  thetas[1] = -M_PI * 0.25;
  //thetas[0] = 0;
  //thetas[1] = 0;

  VECTOR3 edgeBar0 = xBar1 - xBar0;
  VECTOR3 edgeBar1 = xBar2 - xBar1;

  const VECTOR3 x0 = randomVector3();  
  const VECTOR3 x1 = randomVector3();  
  const VECTOR3 x2 = randomVector3();  
  //const VECTOR3 x0 = VECTOR3(0,0,0);
  //const VECTOR3 x1 = VECTOR3(1,0,0);
  //const VECTOR3 x2 = VECTOR3(2,0,0);
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;

  MATRIX3x2 Ebar;
  Ebar.col(0) = edgeBar0;
  Ebar.col(1) = edgeBar1;

  MATRIX3x2 E;
  E.col(0) = edge0;
  E.col(1) = edge1;
  
  const VECTOR3 kbBar = bending.binormal(edgeBar0, edgeBar1);
  const VECTOR3 kb    = bending.binormal(edge0, edge1);

  const VECTOR3 d1_0(0,1,0);
  const VECTOR3 d2_0(0,0,1);
  const VECTOR3 d1_1(0,1,0);
  const VECTOR3 d2_1(0,0,1);

  const VECTOR3 m1_0 =  cos(thetas[0]) * d1_0 + sin(thetas[0]) * d2_0;
  const VECTOR3 m2_0 = -sin(thetas[0]) * d1_0 + cos(thetas[0]) * d2_0;
  const VECTOR3 m1_1 =  cos(thetas[1]) * d1_0 + sin(thetas[1]) * d2_0;
  const VECTOR3 m2_1 = -sin(thetas[1]) * d1_0 + cos(thetas[1]) * d2_0;
  
  MATRIX3x2 M;
  M.col(0) =  0.5 * (m2_0 + m2_1);
  M.col(1) = -0.5 * (m1_0 + m1_1);

  const VECTOR2 kappaBar = M.transpose() * kbBar;
  const VECTOR2 kappa = M.transpose() * kb;

  MATRIX2 B;
  B.setZero();
  B(0,0) = 2.0;
  B(1,1) = 1.0;

  const MATRIX3x2 g0 = bending.edgeGradient(E, B, M, kappa, kappaBar);
  const MATRIX6 H = bending.edgeHessian(E, B, M, kappaBar);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  MATRIX6 finiteDiff;
  while (eps > 1e-8)
  {
    finiteDiff.setZero();

    // for each of the degrees of the freedom
    int index = 0;
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 3; i++, index++)
      {
        // add the pertubation
        MATRIX3x2 Enew = E;
        Enew(i,j) += eps;
  
        const VECTOR3 kbNew = bending.binormal(Enew);
        const VECTOR2 kappaNew = M.transpose() * kbNew;
        
        // get the new psi
        const MATRIX3x2 g = bending.edgeGradient(Enew, B, M, kappaNew, kappaBar);

        // store the finite difference
        finiteDiff.col(index) = (flatten(g) - flatten(g0)) / eps;
      }

    MATRIX6 diff = H - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / H.norm())) / 36.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 36.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
    //cout << " analytic: " << endl << g << endl;
    //cout << " numerical: " << endl << finiteDiff << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the spatial gradient of sin-based bending
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestGradientSinBending(const HOBAK::STRAND::SIN_BENDING& bending)
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand gradient for " << bending.name().c_str() << " bending" << endl;
  cout << "=============================================================== " << endl;

  VECTOR3 xBar0(0,0,0);
  VECTOR3 xBar1(1,1,0);
  VECTOR3 xBar2(2,0,1);

  VECTOR9 positionsBar;
  positionsBar.block<3,1>(0,0) = xBar0;
  positionsBar.block<3,1>(3,0) = xBar1;
  positionsBar.block<3,1>(6,0) = xBar2;

  VECTOR2 thetas;
  thetas[0] =  M_PI * 0.25;
  thetas[1] = -M_PI * 0.25;

  VECTOR3 edgeBar0 = xBar1 - xBar0;
  VECTOR3 edgeBar1 = xBar2 - xBar1;
 
  // TODO: noise this up
  const VECTOR3 x0 = randomVector3();  
  const VECTOR3 x1 = randomVector3();  
  const VECTOR3 x2 = randomVector3();  
  //const VECTOR3 x0 = xBar0 + VECTOR3(0,1,0);
  //const VECTOR3 x1 = xBar1;
  //const VECTOR3 x2 = xBar2;
  VECTOR9 positions;
  positions.block<3,1>(0,0) = x0;
  positions.block<3,1>(3,0) = x1;
  positions.block<3,1>(6,0) = x2;
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;

  MATRIX3x2 Ebar;
  Ebar.col(0) = edgeBar0;
  Ebar.col(1) = edgeBar1;

  MATRIX3x2 E;
  E.col(0) = edge0;
  E.col(1) = edge1;
  
  const VECTOR3 kbBar = bending.binormal(edgeBar0, edgeBar1);
  const VECTOR3 kb    = bending.binormal(edge0, edge1);

  const VECTOR3 d1_0(0,1,0);
  const VECTOR3 d2_0(0,0,1);
  const VECTOR3 d1_1(0,1,0);
  const VECTOR3 d2_1(0,0,1);

  const VECTOR3 m1_0 =  cos(thetas[0]) * d1_0 + sin(thetas[0]) * d2_0;
  const VECTOR3 m2_0 = -sin(thetas[0]) * d1_0 + cos(thetas[0]) * d2_0;
  const VECTOR3 m1_1 =  cos(thetas[1]) * d1_0 + sin(thetas[1]) * d2_0;
  const VECTOR3 m2_1 = -sin(thetas[1]) * d1_0 + cos(thetas[1]) * d2_0;
  
  MATRIX3x2 M;
  M.col(0) =  0.5 * (m2_0 + m2_1);
  M.col(1) = -0.5 * (m1_0 + m1_1);

  const VECTOR2 kappaBar = M.transpose() * kbBar;
  const VECTOR2 kappa = M.transpose() * kb;

  MATRIX2 B;
  B.setZero();
  B(0,0) = 2.0;
  B(1,1) = 1.0;

  MATRIX3x2 M0;
  M0.col(0) = m1_0;
  M0.col(1) = m2_0;
    
  MATRIX3x2 M1;
  M1.col(0) = m1_1;
  M1.col(1) = m2_1;

  const REAL psi0 = bending.psi(B, kappa, kappaBar);
  const VECTOR11 g = bending.gradient(E, B, M0, M1, M, kappa, kappaBar);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  VECTOR11 finiteDiff;
  while (eps > 1e-8)
  {
    finiteDiff.setZero();

    // position degrees
    for (int i = 0; i < 9; i++)
    {
      VECTOR9 positionsNew = positions;
      positionsNew[i] += eps;
  
      const VECTOR3 x0New = positionsNew.segment<3>(0);
      const VECTOR3 x1New = positionsNew.segment<3>(3);
      const VECTOR3 x2New = positionsNew.segment<3>(6);

      const VECTOR3 edge0New = x1New - x0New;
      const VECTOR3 edge1New = x2New - x1New;
      MATRIX3x2 Enew;
      Enew.col(0) = edge0New;
      Enew.col(1) = edge1New;

      const VECTOR3 kbNew = bending.binormal(Enew);
      const VECTOR2 kappaNew = M.transpose() * kbNew;
      
      // get the new psi
      double psi = bending.psi(B, kappaNew, kappaBar);

      // store the finite difference
      finiteDiff[i] = (psi - psi0) / eps;
    }
    for (int i = 0; i < 2; i++)
    {
      // add the pertubation
      VECTOR2 thetasNew = thetas;
      thetasNew[i] += eps;

      const VECTOR3 m1_0new =  cos(thetasNew[0]) * d1_0 + sin(thetasNew[0]) * d2_0;
      const VECTOR3 m2_0new = -sin(thetasNew[0]) * d1_0 + cos(thetasNew[0]) * d2_0;
      const VECTOR3 m1_1new =  cos(thetasNew[1]) * d1_0 + sin(thetasNew[1]) * d2_0;
      const VECTOR3 m2_1new = -sin(thetasNew[1]) * d1_0 + cos(thetasNew[1]) * d2_0;
      MATRIX3x2 Mnew;
      Mnew.col(0) =  0.5 * (m2_0new + m2_1new);
      Mnew.col(1) = -0.5 * (m1_0new + m1_1new);
      const VECTOR2 kappaNew = Mnew.transpose() * kb;
      
      // get the new psi
      double psi = bending.psi(B, kappaNew, kappaBar);

      // store the finite difference
      finiteDiff[9 + i] = (psi - psi0) / eps;
    }

    VECTOR11 thetasLast = finiteDiff;

    finiteDiff.segment<3>(0) = thetasLast.segment<3>(0);
    finiteDiff.segment<3>(4) = thetasLast.segment<3>(3);
    finiteDiff.segment<3>(8) = thetasLast.segment<3>(6);
    finiteDiff[3] = thetasLast[9];
    finiteDiff[7] = thetasLast[10];

    VECTOR11 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 11.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 11.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " g: " << endl << g << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
    //cout << " analytic: " << endl << g << endl;
    //cout << " numerical: " << endl << finiteDiff << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the spatial gradient of Discrete Viscous Threads
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestGradientDVT()
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand spatial gradient for Discrete Viscous Threads" << endl;
  cout << "=============================================================== " << endl;

  BERGOU_2010 DVT(1.0);

  VECTOR9 positions;
  VECTOR3 x0(0,0,0);
  VECTOR3 x1(1,0,0);
  VECTOR3 x2(1,1,0);
  //VECTOR3 x0 = randomVector3(2.0);
  //VECTOR3 x1 = randomVector3(2.0);
  //VECTOR3 x2 = randomVector3(2.0);

  positions.block<3,1>(0,0) = x0;
  positions.block<3,1>(3,0) = x1;
  positions.block<3,1>(6,0) = x2;

  //VECTOR3 edge0 = x0 - x1;
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;
  
  //VECTOR3 tangent0 = edge0 / edge0.norm();
  //VECTOR3 tangent1 = edge1 / edge1.norm();

  const VECTOR2 kappaBar(0,-2);
  const VECTOR2 kappa(1,-1);

  const VECTOR3 director1_0(0,1,0);
  const VECTOR3 director1_1(0,0,1);
  const VECTOR3 director2_0(0,0,1);
  const VECTOR3 director2_1(1,0,0);

  const REAL len = 1.0;
  const VECTOR3 kb(0,0,2);

  BERGOU_2010::_edge0 = edge0;
  BERGOU_2010::_edge1 = edge1;

  /*
  {
  const VECTOR3 sum1 = director1_0 + director1_1;
  const VECTOR3 sum2 = director2_0 + director2_1;
  cout << " Kappa guess: " << 0.5 * VECTOR2(kb.dot(sum2), -kb.dot(sum1)) << endl;
  }
  */

  BERGOU_2010::_kappaBar = kappaBar;
  BERGOU_2010::_kappa    = kappa;
  BERGOU_2010::_len      = len;
  BERGOU_2010::_B        = MATRIX2::Identity();
  BERGOU_2010::_director1_0 = director1_0;
  BERGOU_2010::_director2_0 = director2_0;
  BERGOU_2010::_director1_1 = director1_1;
  BERGOU_2010::_director2_1 = director2_1;
  BERGOU_2010::_kb          = kb;

  REAL psi0 = DVT.psi();
  VECTOR9 g = DVT.gradient();

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR9 finiteDiff;

    // for each of the degrees of the freedom
    int index = 0;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++, index++)
      {
        // add the pertubation
        VECTOR9 positionsNew = positions;
        positionsNew[index] += eps;

        // retrieve the positions
        x0 = positionsNew.block<3,1>(0,0);
        x1 = positionsNew.block<3,1>(3,0);
        x2 = positionsNew.block<3,1>(6,0);
        //VECTOR3 edge0new = x0 - x1;
        VECTOR3 edge0new = x1 - x0;
        VECTOR3 edge1new = x2 - x1;
   
        VECTOR3 tangent0new = edge0new / edge0new.norm();
        VECTOR3 tangent1new = edge1new / edge1new.norm();

        //VECTOR3 director1_0new = parallelTransport(tangent0, tangent0new, director1_0);
        //VECTOR3 director2_0new = parallelTransport(tangent0, tangent0new, director2_0);
        //VECTOR3 director1_1new = parallelTransport(tangent1, tangent1new, director1_1);
        //VECTOR3 director2_1new = parallelTransport(tangent1, tangent1new, director2_1);
        VECTOR3 director1_0new = director1_0;
        VECTOR3 director2_0new = director2_0;
        VECTOR3 director1_1new = director1_1;
        VECTOR3 director2_1new = director2_1;

        const VECTOR3 kb = curvatureBinormal(tangent0new, tangent1new);
        const VECTOR3 sum1 = director1_0new + director1_1new;
        const VECTOR3 sum2 = director2_0new + director2_1new;
        //BERGOU_2010::_director1_0 = director1_0new;
        //BERGOU_2010::_director2_0 = director2_0new;
        //BERGOU_2010::_director1_1 = director1_1new;
        //BERGOU_2010::_director2_1 = director2_1new;
        BERGOU_2010::_kb    = kb;
        BERGOU_2010::_kappa = 0.5 * VECTOR2(kb.dot(sum2), -kb.dot(sum1));
        //cout << " kappa: " << (BERGOU_2010::_kappa).transpose() << endl;

        // get the new edge information
        BERGOU_2010::_edge0 = edge0new;
        BERGOU_2010::_edge1 = edge1new;

        // get the new psi
        double psi = DVT.psi();

        // store the finite difference
        finiteDiff[index] = (psi - psi0) / eps;
      }

    VECTOR9 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 9.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 9.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " g: " << endl << g.transpose() << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff).transpose() << endl;
      cout << " diff: " << endl << clampSmalls(diff).transpose() << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the spatial hessian of Discrete Viscous Threads
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessianDVT()
{
  cout.flush();
  using namespace HOBAK::STRAND;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand spatial hessian for Discrete Viscous Threads" << endl;
  cout << "=============================================================== " << endl;

  BERGOU_2010 DVT(1.0);

  VECTOR9 positions;
  VECTOR3 x0(0,0,0);
  VECTOR3 x1(1,0,0);
  VECTOR3 x2(1,1,0);

  positions.block<3,1>(0,0) = x0;
  positions.block<3,1>(3,0) = x1;
  positions.block<3,1>(6,0) = x2;
  cout << "positions: " << positions.transpose() << endl;

  //VECTOR3 edge0 = x0 - x1;
  VECTOR3 edge0 = x1 - x0;
  VECTOR3 edge1 = x2 - x1;
  
  //VECTOR3 tangent0 = edge0 / edge0.norm();
  //VECTOR3 tangent1 = edge1 / edge1.norm();

  const VECTOR2 kappaBar(0,-2);
  const VECTOR2 kappa(1,-1);

  const VECTOR3 director1_0(0,1,0);
  const VECTOR3 director1_1(0,0,1);
  const VECTOR3 director2_0(0,0,1);
  const VECTOR3 director2_1(1,0,0);

  const REAL len = 1.0;
  const VECTOR3 kb(0,0,2);

  BERGOU_2010::_edge0 = edge0;
  BERGOU_2010::_edge1 = edge1;

  BERGOU_2010::_kappaBar = kappaBar;
  BERGOU_2010::_kappa    = kappa;
  BERGOU_2010::_len      = len;
  BERGOU_2010::_B        = MATRIX2::Identity();
  BERGOU_2010::_director1_0 = director1_0;
  BERGOU_2010::_director2_0 = director2_0;
  BERGOU_2010::_director1_1 = director1_1;
  BERGOU_2010::_director2_1 = director2_1;
  BERGOU_2010::_kb          = kb;
  const VECTOR3 sum1 = director1_0 + director1_1;
  const VECTOR3 sum2 = director2_0 + director2_1;

  VECTOR9 g0 = DVT.gradient();
  MATRIX9 H = DVT.hessian();

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX9 finiteDiff;
    cout << " g0: " << g0.transpose() << endl;

    // for each of the degrees of the freedom
    int index = 0;
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++, index++)
      {
        // add the pertubation
        VECTOR9 positionsNew = positions;
        positionsNew[index] += eps;

        // retrieve the positions
        x0 = positionsNew.block<3,1>(0,0);
        x1 = positionsNew.block<3,1>(3,0);
        x2 = positionsNew.block<3,1>(6,0);
        //cout << " x0: " << x0.transpose() << endl;
        //cout << " x1: " << x1.transpose() << endl;
        //cout << " x2: " << x2.transpose() << endl;
        //VECTOR3 edge0new = x0 - x1;
        VECTOR3 edge0new = x1 - x0;
        VECTOR3 edge1new = x2 - x1;
   
        VECTOR3 tangent0new = edge0new / edge0new.norm();
        VECTOR3 tangent1new = edge1new / edge1new.norm();

        //VECTOR3 director1_0new = parallelTransport(tangent0, tangent0new, director1_0);
        //VECTOR3 director2_0new = parallelTransport(tangent0, tangent0new, director2_0);
        //VECTOR3 director1_1new = parallelTransport(tangent1, tangent1new, director1_1);
        //VECTOR3 director2_1new = parallelTransport(tangent1, tangent1new, director2_1);
        //VECTOR3 director1_0new = director1_0;
        //VECTOR3 director2_0new = director2_0;
        //VECTOR3 director1_1new = director1_1;
        //VECTOR3 director2_1new = director2_1;

        const VECTOR3 kb = curvatureBinormal(tangent0new, tangent1new);
        //const VECTOR3 sum1 = director1_0new + director1_1new;
        //const VECTOR3 sum2 = director2_0new + director2_1new;
        //BERGOU_2010::_director1_0 = director1_0new;
        //BERGOU_2010::_director2_0 = director2_0new;
        //BERGOU_2010::_director1_1 = director1_1new;
        //BERGOU_2010::_director2_1 = director2_1new;
        BERGOU_2010::_kb    = kb;
        BERGOU_2010::_kappa = 0.5 * VECTOR2(kb.dot(sum2), -kb.dot(sum1));
        //cout << " kappa: " << (BERGOU_2010::_kappa).transpose() << endl;

        // get the new edge information
        BERGOU_2010::_edge0 = edge0new;
        BERGOU_2010::_edge1 = edge1new;

        // get the new psi
        VECTOR9 g = DVT.gradient();

        // store the finite difference
        finiteDiff.col(index) = (g - g0) / eps;
      }

    MATRIX9 diff = H - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / H.norm())) / 81.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 81.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H.transpose() << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff).transpose() << endl;
      cout << " diff: " << endl << clampSmalls(diff).transpose() << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the spatial gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestGradient(const HOBAK::STRAND::STRETCHING* material, const vector<VECTOR3>& p)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING strand spatial gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  const REAL dmInv = 1.2345;
  VECTOR3 f = (p[0] - p[1]) * dmInv;

  REAL psi0 = material->psi(f);
  VECTOR6 g = material->spatialGradient(p, dmInv);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR6 finiteDiff;

    // for each of the degrees of the freedom
    int index = 0;
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 3; i++, index++)
      {
        vector<VECTOR3> pNew = p;
        pNew[j][i] += eps;

        // get the new psi
        VECTOR3 fNew = (pNew[0] - pNew[1]) * dmInv;
        double psi = material->psi(fNew);

        // store the finite difference
        finiteDiff[index] = (psi - psi0) / eps;
      }

    VECTOR6 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 6.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 6.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " g: " << endl << g << endl;
      cout << " finite diff: " << endl << finiteDiff << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;

    e++;
  }
  if (relativeMinSeen < 1e-6 || absMinSeen <= 1e-8)
  {
    cout << " TEST PASSED. " << endl;
  }
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the strand bending hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const MATRIX3x2& E)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 PK1 = material->PK1(E, theta0);
  MATRIX6 H = material->hessian(E, theta0);

  const REAL PK1norm = (PK1.norm() > 0.0) ? PK1.norm() : 1.0;

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX6 finiteDiffH;

    // for each of the degrees of the freedom
    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++, i++)
      {
        MATRIX3x2 fnew = E;
        fnew(x,y) += eps;

        // get the new psi
        MATRIX3x2 PK1new = material->PK1(fnew, theta0);

        // store the finite difference
        MATRIX3x2 fd = (PK1new - PK1) / eps;
        finiteDiffH.col(i) = flatten(fd);
      }

    MATRIX6 diff = H - finiteDiffH;

    // cout << " diff: " << endl << diff << endl;
    // cout << " PK1: " << endl << PK1 << endl;
    REAL diffNorm = (fabs(diff.norm() / PK1norm)) / 36.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the strand twisting hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTwistHessian(const HOBAK::VOLUME::TET_STRAND_TWIST* material, const MATRIX3& F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Twist Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  MATRIX9 analyticHess = material->hessian(F);
  MATRIX3 PK1 = material->PK1(F);
  VECTOR9 PK1v = flatten(PK1);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX9 finiteDiff;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++, entry++)
      {
        MATRIX3 vNew = F;
        vNew(j,i) += eps;

        // get the new psi
        MATRIX3 PK1n = material->PK1(vNew);
        VECTOR9 PK1nv = flatten(PK1n);

        // store the finite difference
        finiteDiff.col(entry) = (PK1nv - PK1v) / eps;
      }

    MATRIX9 diff = analyticHess - finiteDiff;
    REAL diffNorm = (abs(diff.norm() / analyticHess.norm())) / 81;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " 2nd F deriv: " << endl << analyticHess << endl;
      cout << " finite diff: " << endl << finiteDiff<< endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}


//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the strand bending hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessianInverted(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const MATRIX3x2& E)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending Inverted Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 PK1 = material->PK1(E, theta0, true);
  MATRIX6 H = material->hessian(E, theta0, true);

  const REAL PK1norm = (PK1.norm() > 0.0) ? PK1.norm() : 1.0;

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX6 finiteDiffH;

    // for each of the degrees of the freedom
    int i = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++, i++)
      {
        MATRIX3x2 fnew = E;
        fnew(x,y) += eps;

        // get the new psi
        MATRIX3x2 PK1new = material->PK1(fnew, theta0, true);

        // store the finite difference
        MATRIX3x2 fd = (PK1new - PK1) / eps;
        finiteDiffH.col(i) = flatten(fd);
      }

    MATRIX6 diff = H - finiteDiffH;

    //cout << " diff: " << endl << diff << endl;
    //cout << " PK1: " << endl << PK1 << endl;
    REAL diffNorm = (fabs(diff.norm() / PK1norm)) / 36.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << H << endl;
      cout << " finite diff: " << endl << finiteDiffH << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// helper for computing a 9x12 dFdX (assuming constant projection coeffs)
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 computeTwistPFpxConst(const REAL& alph, const REAL& bet)
{
  MATRIX9x12 pFpX;
  pFpX.setZero();
  pFpX(0,0) = 1;
  pFpX(0,3) = -1 + alph;
  pFpX(0,6) = -alph;
  pFpX(1,1) = 1;
  pFpX(1,4) = -1 + alph;
  pFpX(1,7) = -alph;
  pFpX(2,2) = 1;
  pFpX(2,5) = -1 + alph;
  pFpX(2,8) = -alph;
  pFpX(3,3) = bet;
  pFpX(3,6) = -1 - bet;
  pFpX(3,9) = 1;
  pFpX(4,4) = bet;
  pFpX(4,7) = -1 - bet;
  pFpX(4,10) = 1;
  pFpX(5,5) = bet;
  pFpX(5,8) = -1 - bet;
  pFpX(5,11) = 1;
  pFpX(6,3) = -1;
  pFpX(6,6) = 1;
  pFpX(7,4) = -1;
  pFpX(7,7) = 1;
  pFpX(8,5) = -1;
  pFpX(8,8) = 1;
  
  return pFpX;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on position-based first derivative (Tet twisting)
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTwistForcesConst(const HOBAK::VOLUME::TET_STRAND_TWIST& material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING TET_STRAND_TWIST forces " << endl;
  cout << "=============================================================== " << endl;
  vector<VECTOR3> vertices;
  vertices.resize(4);

  // randomize the vertices
  for (unsigned int x = 0; x < vertices.size(); x++)
    vertices[x] += randomVector3(1.0);
  
  // cache the current state
  const vector<VECTOR3> v0 = vertices;

  // get F from the vertices
  VECTOR3 e0 = vertices[0] - vertices[1];
  VECTOR3 e1 = vertices[2] - vertices[1];
  VECTOR3 e2 = vertices[3] - vertices[2];
  REAL alph = e0.dot(e1)/e1.squaredNorm();
  e0 = e0 - alph * e1;
  REAL bet = e2.dot(e1)/e1.squaredNorm();
  e2 = e2 - bet * e1;
  MATRIX3 F;
  F.col(0) = e0;
  F.col(1) = e2;
  F.col(2) = e1;

  // get dF/dx (treating it as constant alpha/beta)
  MATRIX9x12 dFdX = computeTwistPFpxConst(alph,bet);

  // get the force
  VECTOR12 forces = -dFdX.transpose() * flatten(material.PK1(F));

  // get the current energy
  REAL psi0 = material.psi(F);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-9)
  {
    VECTOR finiteDiff(forces.size());
    finiteDiff.setZero();

    int index = 0;
    for (unsigned int y = 0; y < vertices.size(); y++)
    {
      for (unsigned int x = 0; x < 3; x++, index++)
      {
        // reset the other components
        vertices[y] = v0[y];

        // just perturb one component
        vertices[y][x] += eps;
       
        // compute the energy
        // get F from the vertices
        VECTOR3 e0 = vertices[0] - vertices[1];
        VECTOR3 e1 = vertices[2] - vertices[1];
        VECTOR3 e2 = vertices[3] - vertices[2];
        REAL alph = e0.dot(e1)/e1.squaredNorm();
        e0 = e0 - alph * e1;
        REAL bet = e2.dot(e1)/e1.squaredNorm();
        e2 = e2 - bet * e1;
        MATRIX3 F;
        F.col(0) = e0;
        F.col(1) = e2;
        F.col(2) = e1;

        REAL psi1 = material.psi(F);

        finiteDiff[index] = -(psi1 - psi0) / eps;
      }

      // last reset
      vertices[y] = v0[y];
    }

    VECTOR diff = forces - finiteDiff;
    REAL diffNorm = fabs(diff.norm() / forces.norm());
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 5 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " forces: " << endl << forces.transpose() << endl;
      cout << " finite diff: " << endl << finiteDiff.transpose() << endl;
      cout << " diff: " << endl << diff.transpose() << endl;
      return false;
    }
    else
      
    eps *= 0.1;
    e++;
  }
  cout << " TEST PASSED. " << endl;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1 (Tet twisting)
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1(const HOBAK::VOLUME::TET_STRAND_TWIST* material, const MATRIX3& F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Twist PK1 for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  MATRIX3 analyticPK1 = material->PK1(F);
  REAL psi0 = material->psi(F);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-10)
  {
    MATRIX3 finiteDiffPK1;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++, entry++)
      {
        MATRIX3 vNew = F;
        vNew(i,j) += eps;

        // get the new psi
        double psi = material->psi(vNew);

        // store the finite difference
        finiteDiffPK1(i,j) = (psi - psi0) / eps;
      }

    MATRIX3 diff = analyticPK1 - finiteDiffPK1;
    REAL diffNorm = (abs(diff.norm() / analyticPK1.norm())) / 9;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 6 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " PK1: " << endl << analyticPK1 << endl;
      cout << " finite diff: " << endl << finiteDiffPK1 << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void computeAlphaBeta(const vector<VECTOR3>& vertices, REAL& alpha, REAL& beta)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
    
  alpha = e0.dot(e1)/e1.squaredNorm();
  beta = e2.dot(e1)/e1.squaredNorm();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX3 computeFtwist(const vector<VECTOR3>& vertices)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  VECTOR3 e0 = vertices[0] - vertices[1];
  VECTOR3 e2 = vertices[3] - vertices[2];
  REAL alpha, beta;
  computeAlphaBeta(vertices, alpha, beta);

  e0 = e0 - alpha * e1;
  e2 = e2 - beta * e1;
  MATRIX3 Ftwist;
  Ftwist.col(0) = e0;
  Ftwist.col(1) = e2;
  Ftwist.col(2) = e1;

  return Ftwist;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 computeTwistPFPx(const vector<VECTOR3>& vertices)
{
  REAL alpha, beta;
  computeAlphaBeta(vertices, alpha, beta);

  MATRIX9x12 pFpX;
  pFpX.setZero();
  pFpX(0,0) = 1;
  pFpX(0,3) = -1 + alpha;
  pFpX(0,6) = -alpha;
  pFpX(1,1) = 1;
  pFpX(1,4) = -1 + alpha;
  pFpX(1,7) = -alpha;
  pFpX(2,2) = 1;
  pFpX(2,5) = -1 + alpha;
  pFpX(2,8) = -alpha;
  pFpX(3,3) = beta;
  pFpX(3,6) = -1 - beta;
  pFpX(3,9) = 1;
  pFpX(4,4) = beta;
  pFpX(4,7) = -1 - beta;
  pFpX(4,10) = 1;
  pFpX(5,5) = beta;
  pFpX(5,8) = -1 - beta;
  pFpX(5,11) = 1;
  pFpX(6,3) = -1;
  pFpX(6,6) = 1;
  pFpX(7,4) = -1;
  pFpX(7,7) = 1;
  pFpX(8,5) = -1;
  pFpX(8,8) = 1;

  return pFpX;
}

//////////////////////////////////////////////////////////////////////////////
// unsimplified original PFPx
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 computeTwistPFPxOrig(const vector<VECTOR3>& vertices)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 e0perp = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 e2perp = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const MATRIX3 z33 = MATRIX3::Zero(), I3 = MATRIX3::Identity();
  const MATRIX3 sigma = I3 - t1*t1.transpose();
  const MATRIX3 eta0 = -(t1.dot(e0)*sigma+t1*e0perp.transpose())/norm1;
  const MATRIX3 eta2 = -(t1.dot(e2)*sigma+t1*e2perp.transpose())/norm1;

  MATRIX9x12 pFpX;
  pFpX << sigma, -eta0-sigma, eta0, z33,
          z33, -eta2, eta2 - sigma, sigma,
          z33, -I3, I3, z33;

  return pFpX;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX12 computeExtraTerm(const vector<VECTOR3>& vertices, const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 z0 = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 z1 = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const VECTOR3 tb = z0.cross(z1).normalized();
  const VECTOR3 tau0 = z0 / normz0, tau1 = z1 / normz1;
  const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized(); 
  const MATRIX3  z33 = MATRIX3::Zero(); const VECTOR3 z3 = VECTOR3::Zero();
  // blocks
  const MATRIX3 t1tau0perp = t1 * tau0perp.transpose(), t1tau1perp = t1 * tau1perp.transpose();
  MATRIX3x12 tau0pSigmapx, tau1pSigmapx, eta2px0, eta3px0, tau0peta2px, tau1peta3px;
  VECTOR12 eta2px1, eta3px1;
  tau0pSigmapx << z33, t1tau0perp, -t1tau0perp, z33;
  tau1pSigmapx << z33, t1tau1perp, -t1tau1perp, z33;
  tau0pSigmapx *= 1.0/norm1; tau1pSigmapx *= 1.0/norm1;
  const VECTOR3 reflect2 = e0 - 2.0 * t1 * (t1.dot(e0)), reflect3 = e2 - 2.0 * t1 * (t1.dot(e2)); 
  eta2px0 << z33, reflect2*tau0perp.transpose(), -reflect2*tau0perp.transpose(), z33;
  eta2px1 << e1, -reflect2 - e1, reflect2, z3;
  eta3px0 << z33, reflect3*tau1perp.transpose(), -reflect3*tau1perp.transpose(), z33;
  eta3px1 << z3, -reflect3, reflect3 - e1, e1;
  tau0peta2px = 1.0/norm1/norm1 *(eta2px0 - tau0perp*eta2px1.transpose());
  tau1peta3px = 1.0/norm1/norm1 *(eta3px0 - tau1perp*eta3px1.transpose());
  // this computes p Psi p theta. should be computed inside the material class.
  const REAL s = tb.dot(t1) > 0?-1:1;
  const REAL thetaGrad = 2.0 * s * material->mu() * (acos(tau0.dot(tau1)) - material->theta0());
  MATRIX12 extra;
  extra << 1.0/normz0 * tau0pSigmapx,
           1.0/normz0*(-tau0peta2px-tau0pSigmapx) + 1.0/normz1 * tau1peta3px,
           1.0/normz0*tau0peta2px - 1.0/normz1*(tau1peta3px - tau1pSigmapx),
          -1.0/normz1*tau1pSigmapx;
  return thetaGrad * extra;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX12 computeExtraTermClamped(const vector<VECTOR3>& vertices, const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  const VECTOR3 e1 = vertices[2] - vertices[1];
  const VECTOR3 e0 = vertices[0] - vertices[1];
  const VECTOR3 e2 = vertices[3] - vertices[2];
  const VECTOR3 z0 = e0 - e0.dot(e1)/e1.squaredNorm()*e1;
  const VECTOR3 z1 = e2 - e2.dot(e1)/e1.squaredNorm()*e1;
  const REAL norm1 = e1.norm();
  const VECTOR3 t1 = e1/norm1;
  const REAL normz0 = z0.norm(), normz1 = z1.norm();
  const VECTOR3 tb = z0.cross(z1).normalized();
  const VECTOR3 tau0 = z0 / normz0, tau1 = z1 / normz1;
  //const VECTOR3 tau0perp = tau0.cross(tb).normalized(); 
  const VECTOR3 tau1perp = tau1.cross(tb).normalized(); 
  //const MATRIX3  z33 = MATRIX3::Zero(); 
  const VECTOR3 z3 = VECTOR3::Zero();

  const REAL alpha = -4.0 * (tau0.dot(tau1perp)) / norm1 / norm1;
  const REAL beta = 2.0 * (tau0.dot(tau1perp))/norm1; 
  const REAL beta2 = beta*beta, alpha2 = alpha*alpha;
  const REAL a1 = t1.dot(e2) / norm1/normz1, b1 = 1/normz1;
  const REAL a0 = t1.dot(e0) / norm1/normz0, b0 = 1/normz0;
  VECTOR12 p0, p1, v0, v1, u0, u1;
  p0 << -b0*t1, (-a1-a0+b0)*t1, (a1+a0+b1)*t1, -b1*t1;
  p1 << b0*t1, (-a1+a0-b0)*t1, (a1-a0+b1)*t1, -b1*t1;
  v0 << z3, tau0, -tau0, z3; v1 << z3, tau1, -tau1, z3; 
  u0 = v0 + v1; u1 = v0 - v1;
  const REAL d00 = (p0.dot(p0)) / (u0.dot(u0)), d01 = (p0.dot(p1)) / (u1.dot(u1));
  const REAL d10 = (p0.dot(p1)) / (u0.dot(u0)), d11 = (p1.dot(p1)) / (u1.dot(u1));
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  SolveP4De(x, -alpha2 - beta2*(d11+d00), beta2*alpha*(d00-d11), beta2*beta2*(d11*d00 - d10*d01));
  MATRIX12 clamped; clamped.setZero();
  // this computes p Psi p theta. should be computed inside the material class.
  const REAL s = tb.dot(t1) > 0?-1:1;
  const REAL thetaGrad = 2.0 * s * material->mu() * (acos(tau0.dot(tau1)) - material->theta0());
  for(int i = 0; i<4; i++){
    if(x[i] * thetaGrad > 0.0){
      const REAL f = beta/x[i], e = (1.0 /f/f - d00 + alpha/beta/f)/d10;
      const VECTOR12 eigVec = u0 * 1.0 + u1 * e + p0 * f + p1 * (e*f);
      clamped += thetaGrad * x[i] * eigVec * eigVec.transpose(); 
    }
  }
  return clamped;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9x12 computeDelTwistPFPx(const vector<VECTOR3>& vertices)
{
  REAL alpha,beta;
  computeAlphaBeta(vertices, alpha, beta);
  VECTOR3 e0 = vertices[0]-vertices[1];
  VECTOR3 e1 = vertices[2]-vertices[1];
  VECTOR3 e2 = vertices[3]-vertices[2];

  // all the nontrivial del term blocks (with 1/e1^2 scalar factored out)
  MATRIX3 De0px0 = -e1*e1.transpose(); ///also De2px3
  MATRIX3 De0px1 = -e1*((2*alpha - 1)*e1 - e0).transpose();
  MATRIX3 De0px2 = -e1*(e0 - 2*alpha*e1).transpose();
  MATRIX3 De2px1 = -e1*(2*beta*e1 - e2).transpose();
  MATRIX3 De2px2 = e1*((2*beta + 1)*e1 - e2).transpose();

  MATRIX9x12 delpFpX;
  delpFpX.setZero();
  delpFpX.block<3,3>(0,0) = De0px0;
  delpFpX.block<3,3>(0,3) = De0px1;
  delpFpX.block<3,3>(0,6) = De0px2;
  delpFpX.block<3,3>(3,3) = De2px1;
  delpFpX.block<3,3>(3,6) = De2px2;
  delpFpX.block<3,3>(3,9) = De0px0;

  return 1/e1.squaredNorm()*delpFpX;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1 (Tet twisting)
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTwistingSpatialGradient(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Twist gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  vector<VECTOR3> vertices0(4);
  for (int x = 0; x < 4; x++)
    vertices0[x] = randomVector3();
  MATRIX3 F0 = computeFtwist(vertices0);

  MATRIX3 analyticPK1 = material->PK1(F0);
  VECTOR12 forceAnalytic = computeTwistPFPx(vertices0).transpose() * flatten(analyticPK1);
  REAL psi0 = material->psi(F0);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-10)
  {
    VECTOR12 forceNumerical;

    // for each of the degrees of the freedom
    int entry = 0;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++, entry++)
      {
        vector<VECTOR3> vertices1 = vertices0;
        vertices1[i][j] += eps;
        const MATRIX3 F1 = computeFtwist(vertices1);

        // get the new psi
        double psi = material->psi(F1);

        // store the finite difference
        forceNumerical[entry] = (psi - psi0) / eps;
      }

    VECTOR12 diff = forceAnalytic - forceNumerical;;
    REAL diffNorm = (abs(diff.norm() / forceAnalytic.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 6 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " analytic: " << endl << forceAnalytic << endl;
      cout << " numerical: " << endl << forceNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const MATRIX3x2& E)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(E, theta0);
  MATRIX3x2 p = material->PK1(E, theta0);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX3x2 finiteDiffP;

    // for each of the degrees of the freedom
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++)
      {
        MATRIX3x2 fnew = E;
        fnew(x,y) += eps;

        // get the new psi
        double psi = material->psi(fnew, theta0);

        // store the finite difference
        finiteDiffP(x,y) = (psi - psi0) / eps;
      }

    MATRIX3x2 diff = p - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / p.norm())) / 6.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " p: " << endl << p << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1, with inversion
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1Inverted(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const MATRIX3x2& E)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending Inverted P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  bool inverted = true;
  REAL psi0 = material->psi(E, theta0, inverted);
  MATRIX3x2 p = material->PK1(E, theta0, inverted);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX3x2 finiteDiffP;

    // for each of the degrees of the freedom
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++)
      {
        MATRIX3x2 fnew = E;
        fnew(x,y) += eps;

        // get the new psi
        double psi = material->psi(fnew, theta0, inverted);

        // store the finite difference
        finiteDiffP(x,y) = (psi - psi0) / eps;
      }

    MATRIX3x2 diff = p - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / p.norm())) / 6.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " p: " << endl << p << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the PK1
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestPK1(const HOBAK::STRAND::STRETCHING* material, const VECTOR3& f)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(f);
  VECTOR3 p = material->PK1(f);

  cout << " psi0: " << psi0 << endl;

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR3 finiteDiffP;

    // for each of the degrees of the freedom
    for (int x = 0; x < 3; x++)
    {
      VECTOR3 fnew = f;
      fnew[x] += eps;

      // get the new psi
      double psi = material->psi(fnew);

      // store the finite difference
      finiteDiffP[x] = (psi - psi0) / eps;
    }

    VECTOR3 diff = p - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / p.norm())) / 3.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " p: " << endl << p << endl;
      cout << " finite diff: " << endl << finiteDiffP << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testPathological()
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Pathological case for QUADRATIC_BENDING" << endl;
  cout << "=============================================================== " << endl;
  using namespace HOBAK::STRAND;
  QUADRATIC_BENDING quadratic(1.0);
  
  MATRIX3x2 pathological;
  pathological << -1,0,
                   0,1,
                   0,0;

  MATRIX6 H = quadratic.hessian(pathological, theta0);
  MATRIX6 Hnumerical = clampEigenvalues(H);

  MATRIX6 Hanalytic = quadratic.clampedHessian(pathological, theta0);

  MATRIX6 diff = Hnumerical - Hanalytic;

  const REAL diffNorm = diff.norm() / 36.0;

  cout << " Numerical: " << endl;
  cout << Hnumerical << endl;
  cout << " Analytic: " << endl;
  cout << Hanalytic << endl;
  cout << " Diff: " << endl;
  cout << diff << endl;
  cout << " Diff norm: " << diffNorm << endl;

  if (diffNorm < 1e-8)
    return true;

  return false;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testClampedHessian(const HOBAK::STRAND::STRETCHING* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  VECTOR3 f = randomVector3();

  MATRIX3 original, clamped, diff;
  REAL diffNorm;

  // verify that with all positive, we still get the same Hessian
  original = clampEigenvalues(material->hessian(f));
  clamped  = material->clampedHessian(f);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " All positive test: PASSED " << endl;
  else
  {
    cout << " All positive test: FAILED " << endl;
    cout << " diff: " << endl;
    cout << diff << endl;
    cout << " Numerical Clamp: " << endl;
    cout << original << endl;
    cout << " Analytic Clamp: " << endl;
    cout << clamped << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    return false;
  }
  
  // try lots of Fs here, as this case tends to be tricky
  std::mt19937 gen(314159);
  std::uniform_real_distribution<REAL> dist(0.0, 1.0);
  for (int x = 0; x < 20; x++)
  {
    // make sure it's in the compression regime
    f = randomVector3(1.0);
    f = f / f.norm();
    f *= dist(gen);

    original = clampEigenvalues(material->hessian(f));
    clamped  = material->clampedHessian(f);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " Compression test: PASSED " << endl;
    else
    {
      MATRIX3 unclamped = material->hessian(f);
      cout << " Compression test: FAILED " << endl;
      cout << " Unclamped: " << endl << unclamped << endl;
      cout << " Numerical clamped: " << endl << original << endl;
      cout << " Analytic clamped: " << endl << clamped << endl;
      cout << " Diff: " << endl << diff << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      cout << " Unclamped eigs: " << eigenvalues(unclamped).transpose() << endl;
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testClampedHessian(const HOBAK::STRAND::ISOTROPIC_BENDING* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  // try a whole bunch of them
  for (int x = 0; x < 100; x++)
  {
    MATRIX3x2 E = randomMatrix3x2();
    MATRIX6 original, clamped, diff;
    REAL diffNorm;

    // verify that with all positive, we still get the same Hessian
    original = clampEigenvalues(material->hessian(E, theta0));
    clamped  = material->clampedHessian(E, theta0);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " All positive test: PASSED " << endl;
    else
    {
      cout << " All positive test: FAILED " << endl;
      cout << " diff: " << endl;
      cout << diff << endl;
      cout << " Numerical Clamp: " << endl;
      cout << original << endl;
      cout << " Analytic Clamp: " << endl;
      cout << clamped << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      return false;
    }
  }
 
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works (Tet strand twist)
//////////////////////////////////////////////////////////////////////////////
bool testClampedHessian(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  // try a whole bunch of them
  for (int x = 0; x < 100; x++)
  {
    MATRIX3 F;
    VECTOR3 e0 = randomVector3();
    VECTOR3 e1 = randomVector3();
    VECTOR3 e2 = randomVector3();
    VECTOR3 e1Hat = e1.normalized();
    e0 = e0 - e0.dot(e1Hat) * e1Hat;
    e2 = e2 - e2.dot(e1Hat) * e1Hat;
    F.col(0) = e0;
    F.col(1) = e2;
    F.col(2) = e1;

    MATRIX9 original, clamped, diff;
    REAL diffNorm;

    // verify that with all positive, we still get the same Hessian
    original = clampEigenvalues(material->hessian(F));
    clamped  = material->clampedHessian(F);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " All positive test: PASSED " << endl;
    else
    {
      cout << " All positive test: FAILED " << endl;
      cout << " diff: " << endl;
      cout << diff << endl;
      cout << " Numerical Clamp: " << endl;
      cout << original << endl;
      cout << " Analytic Clamp: " << endl;
      cout << clamped << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      return false;
    }
  }
 
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works (Tet strand twist)
//////////////////////////////////////////////////////////////////////////////
bool testClampedForceGradient(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  // try a whole bunch of them
  for (int x = 0; x < 100; x++)
  {
    vector<VECTOR3> vertices(4);
    for (int y = 0; y < 4; y++)
      vertices[y] = randomVector3();

    // brute-force it
    const MATRIX3 F = computeFtwist(vertices);
    const MATRIX9x12 pFpX = computeTwistPFPx(vertices);
    const MATRIX9x12 delpFpX = computeDelTwistPFPx(vertices);
    const MATRIX9 H = material->hessian(F);
    const MATRIX9 plusH = clampEigenvalues(H);
    const MATRIX9 minusH = -1.0 * clampEigenvalues(MATRIX9(-1.0 * H));
    const MATRIX12 unfiltered = pFpX.transpose() * H * pFpX - delpFpX.transpose() * H * delpFpX;
    const MATRIX12 original = pFpX.transpose() * plusH * pFpX - delpFpX.transpose() * minusH * delpFpX;
    //const MATRIX12 original = pFpX.transpose() * plusH * pFpX;

    // verify that with all positive, we still get the same Hessian
    const MATRIX12 analytic  = material->clampedForceGradient(vertices);
    const MATRIX12 diff = original - analytic;
    const REAL diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-6 || original.norm() < 1e-8)
      cout << " All positive test: PASSED " << endl;
    else
    {
      cout << " All positive test: FAILED " << endl;
      cout << " diff: " << endl;
      cout << diff << endl;
      cout << " Numerical Clamp: " << endl;
      cout << original << endl;
      cout << " Analytic Clamp: " << endl;
      cout << analytic << endl;

      cout << " plusH: " << eigenvalues(plusH) << endl;
      cout << " minusH: " << eigenvalues(minusH) << endl;
      cout << " unfiltered:     " << eigenvalues(unfiltered).transpose() << endl;
      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(analytic).transpose() << endl;
      return false;
    }
    cout << " unfiltered:     " << eigenvalues(unfiltered).transpose() << endl;
    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(analytic).transpose() << endl;
  }
 
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works (Tet strand twist)
//////////////////////////////////////////////////////////////////////////////
bool testClampedForceGradientRank4(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << " with rank-4 correction" << endl;
  cout << "=============================================================== " << endl;

  // try a whole bunch of them
  for (int x = 0; x < 100; x++)
  {
    vector<VECTOR3> vertices(4);
    for (int y = 0; y < 4; y++)
      vertices[y] = randomVector3();

    // brute-force it
    const MATRIX3 F = computeFtwist(vertices);
    const MATRIX9x12 pFpX = computeTwistPFPxOrig(vertices);
    const MATRIX9 H = clampEigenvalues(material->hessian(F));
    const MATRIX12 unfiltered = material->forceGradientRank4(vertices);
    const MATRIX12 original = pFpX.transpose() * H * pFpX + clampEigenvalues(computeExtraTerm(vertices, material));

    // verify that with all positive, we still get the same Hessian
    const MATRIX12 analytic  = material->clampedForceGradientRank4(vertices);
    const MATRIX12 diff = original - analytic;
    const REAL diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t" << endl;

    /*
    if (diffNorm < 1e-6 || original.norm() < 1e-8)
      cout << " All positive test: PASSED " << endl;
    else
    {
      cout << " All positive test: FAILED " << endl;
      cout << " diff: " << endl;
      cout << diff << endl;
      cout << " Numerical Clamp: " << endl;
      cout << original << endl;
      cout << " Analytic Clamp: " << endl;
      cout << analytic << endl;

      cout << " unfiltered:     " << clampSmalls(eigenvalues(unfiltered).transpose()) << endl;
      cout << " Numerical eigs: " << clampSmalls(eigenvalues(original).transpose()) << endl;
      cout << " Analytic eigs:  " << clampSmalls(eigenvalues(analytic).transpose()) << endl;
      return false;
    }
    */
    cout << " unfiltered:     " << clampSmalls(eigenvalues(unfiltered).transpose()) << endl;
    cout << " Numerical eigs: " << clampSmalls(eigenvalues(original).transpose()) << endl;
    cout << " Analytic eigs:  " << clampSmalls(eigenvalues(analytic).transpose()) << endl;
  }
 
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testSpatialClampedHessian(const HOBAK::STRAND::STRETCHING* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING spatial eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  vector<VECTOR3> p(2);
  p[0] = randomVector3();
  p[1] = randomVector3();
  const REAL dmInv = 1.12345;

  MATRIX6 original, clamped, diff;
  REAL diffNorm;

  // verify that with all positive, we still get the same Hessian
  original = clampEigenvalues(material->spatialHessian(p, dmInv));
  clamped  = material->spatialClampedHessian(p, dmInv);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " All positive test: PASSED " << endl;
  else
  {
    cout << " All positive test: FAILED " << endl;
    cout << " diff: " << endl;
    cout << diff << endl;
    cout << " Numerical Clamp: " << endl;
    cout << original << endl;
    cout << " Analytic Clamp: " << endl;
    cout << clamped << endl;

    cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
    cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
    return false;
  }
 
  // try lots of Fs here, as this case tends to be tricky
  std::mt19937 gen(314159);
  std::uniform_real_distribution<REAL> dist(0.0, 1.0);
  for (int x = 0; x < 20; x++)
  {
    // make sure it's in the compression regime
    vector<VECTOR3> p(2);
    p[0] = randomVector3();
    p[1] = randomVector3();
    VECTOR3 d = p[1] - p[0];
    d = d / d.norm();
    p[1] = p[0] + dist(gen) * d;

    const REAL unitDmInv = 1.0;

    original = clampEigenvalues(material->spatialHessian(p, unitDmInv));
    clamped  = material->spatialClampedHessian(p, unitDmInv);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " Compression test: PASSED " << endl;
    else
    {
      MATRIX6 unclamped = material->spatialHessian(p, unitDmInv);
      cout << " Compression test: FAILED " << endl;
      cout << " Unclamped: " << endl << unclamped << endl;
      cout << " Numerical clamped: " << endl << original << endl;
      cout << " Analytic clamped: " << endl << clamped << endl;
      cout << " Diff: " << endl << diff << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      cout << " Unclamped eigs: " << eigenvalues(unclamped).transpose() << endl;
      return false;
    }
  }
  return true;
}
#if 0
//////////////////////////////////////////////////////////////////////////////
// simple twisting tests based on known solutions
//////////////////////////////////////////////////////////////////////////////
void testTwisting()
{
  // procedurally generate a scene
  std::vector<VECTOR3> restVertices;
  std::vector<VECTOR3> vertices;
  const int maxElements = 5;
  const REAL frequency = 1.0;
  VECTOR restThetas(maxElements - 1);
  VECTOR thetas(maxElements - 1);
  for (int x = 0; x < maxElements; x++)
  {
    double frac = (double)x / maxElements;
    restVertices.push_back(VECTOR3(cos(frequency * frac * 2 * M_PI), 
                                  sin(frequency * frac * 2 * M_PI), 
                                  3.0 * frac));
    vertices.push_back(VECTOR3(2.0 *cos(frequency * frac * 2 * M_PI), 
                               3.0 * sin(frequency * frac * 2 * M_PI), 
                               4.0 * 3.0 * frac));
    if (x != maxElements - 1)
    {
      restThetas[x] = 0.0;
      thetas[x] = M_PI * frac;
    }
  }

  // initialize the strand mesh
  STRAND_MESH strandMesh(restVertices, restThetas);
  strandMesh.setPositions(vertices);

  STRAND::QUADRATIC_STRETCHING stretchingEnergy(1.0);
  const REAL currentStretch = strandMesh.computeStretchingEnergy(stretchingEnergy);
  cout << " Stretching energy: " << currentStretch << endl;

  cout << " Current positions: " << endl;
  cout << strandMesh.getPositions() << endl;
  cout << " Current twists: " << endl;
  cout << strandMesh.twistAngles() << endl;

  cout.precision(16);
  const REAL currentTwisting = strandMesh.computeTwistingEnergy(vertices, restThetas);
  cout << " Twisting energy after update 0: " << currentTwisting << endl;
  cout << " Should be:                      0.05803780184767759" << endl;
  REQUIRE(fabs(currentTwisting - 0.05803780184767759) < 1e-8);
  cout << "============================================== " << endl;
  cout << " TEST 1 PASSED" << endl;
  cout << "============================================== " << endl;

  vector<VECTOR3> points(5);
  vector<VECTOR3> pointsOld;
  vector<VECTOR3> directorOld;
  pointsOld = vertices;
  directorOld = strandMesh.director1();

  points[0] = VECTOR3(2, 0, 0);
  points[1] = VECTOR3(0.6180339887498949, 2.85316954888546, 2.4);
  points[2] = VECTOR3(-1.618033988749895, 2.580046404216302, 3.002379893500697);
  points[3] = VECTOR3(-1.652647009869095, 1.468123264721223, 3.714122257397702);
  points[4] = VECTOR3(-0.7397848211988234, 1.081774310225489, 4.587970986481197);
  thetas[0] = -0.0378268354684716;
  thetas[1] = -0.04544452057262368;
  thetas[2] = 0.03516919694788014;
  thetas[3] = -0.08891637502264385;

  REAL currentTwisting2 = 0;
  for (unsigned int x = 0; x < points.size() - 2; x++)
  {
    vector<VECTOR3> positions;
    vector<VECTOR3> positionsOld;
    vector<VECTOR3> directorsOld;
    vector<REAL> restEdgeLengths;
    vector<REAL> newThetas;

    positions.push_back(points[x]);
    positions.push_back(points[x + 1]);
    positions.push_back(points[x + 2]);
    positionsOld.push_back(pointsOld[x]);
    positionsOld.push_back(pointsOld[x + 1]);
    positionsOld.push_back(pointsOld[x + 2]);
    directorsOld.push_back(directorOld[x]);
    directorsOld.push_back(directorOld[x + 1]);

    restEdgeLengths.push_back(strandMesh.restEdgeLengths()[x]);
    restEdgeLengths.push_back(strandMesh.restEdgeLengths()[x + 1]);

    REAL restTwist = strandMesh.restTwistAngles()[x];

    newThetas.push_back(thetas[x]);
    newThetas.push_back(thetas[x + 1]);

    currentTwisting2 += strandMesh.computeTwistingEnergy(positions, positionsOld, directorsOld, restEdgeLengths, restTwist, newThetas);
  }
  strandMesh.setPositions(points);
  strandMesh.twistAngles() = thetas;
  strandMesh.updateAngles();
  cout << endl;
  cout << " Twisting energy after update 1: " << currentTwisting2 << endl;
  cout << " Should be:                      0.001720924626346224" << endl;
  REQUIRE(fabs(currentTwisting2 - 0.001720924626346224) < 1e-8);
  cout << "============================================== " << endl;
  cout << " TEST 2 PASSED" << endl;
  cout << "============================================== " << endl;
 
  pointsOld = points;
  directorOld = strandMesh.director1();
  points[0] = VECTOR3(2, 0, 0);
  points[1] = VECTOR3(0.6180339887498949, 2.85316954888546, 2.4);
  points[2] = VECTOR3(-1.618033988749895, 2.601024953596726, 2.945981555886704);
  points[3] = VECTOR3(-1.714556505972326, 1.460882431264167, 3.594458586846398);
  points[4] = VECTOR3(-0.8244267850991962, 1.015924414217568, 4.468376168745999);
  thetas[0] = -0.1015788648851949;
  thetas[1] = -0.08855006995666639;
  thetas[2] = 0.04813479563490996;
  thetas[3] = -0.0609632311622805;

  REAL currentTwisting3 = 0;
  for (unsigned int x = 0; x < points.size() - 2; x++)
  {
    vector<VECTOR3> positions;
    vector<VECTOR3> positionsOld;
    vector<VECTOR3> directorsOld;
    vector<REAL> restEdgeLengths;
    vector<REAL> newThetas;

    positions.push_back(points[x]);
    positions.push_back(points[x + 1]);
    positions.push_back(points[x + 2]);
    positionsOld.push_back(pointsOld[x]);
    positionsOld.push_back(pointsOld[x + 1]);
    positionsOld.push_back(pointsOld[x + 2]);
    directorsOld.push_back(directorOld[x]);
    directorsOld.push_back(directorOld[x + 1]);

    restEdgeLengths.push_back(strandMesh.restEdgeLengths()[x]);
    restEdgeLengths.push_back(strandMesh.restEdgeLengths()[x + 1]);

    REAL restTwist = strandMesh.restTwistAngles()[x];

    newThetas.push_back(thetas[x]);
    newThetas.push_back(thetas[x + 1]);

    currentTwisting3 += strandMesh.computeTwistingEnergy(positions, positionsOld, directorsOld, restEdgeLengths, restTwist, newThetas);
  }
  strandMesh.setPositions(points);
  strandMesh.twistAngles() = thetas;
  strandMesh.updateAngles();

  cout << " Twisting energy after update 3: " << currentTwisting3 << endl;
  cout << " Should be:                      0.001528117491321596" << endl;
  REQUIRE(fabs(currentTwisting3 - 0.001528117491321596) < 1e-8);
  cout << "============================================== " << endl;
  cout << " TEST 3 PASSED" << endl;
  cout << "============================================== " << endl;
}
#endif

//////////////////////////////////////////////////////////////////////////////
// run a small strand regression scene
//////////////////////////////////////////////////////////////////////////////
void strandRegressionScene()
{
  using namespace STRAND;

  // procedurally generate a scene
  std::vector<VECTOR3> restVertices;

  // L shape
  restVertices.push_back(VECTOR3(-1,0,0));
  restVertices.push_back(VECTOR3(0,0,0));
  restVertices.push_back(VECTOR3(1,0,0));
  restVertices.push_back(VECTOR3(1,0,-1));

  vector<VECTOR3> vertices = restVertices;
  vertices[2] = VECTOR3(0,1,0);
  vertices[3] = VECTOR3(1,1,0);

  const REAL E = 10;
  const REAL nu = 0.025 - 1.0;  // non-physical, but good for debugging
  const REAL density = 1e-3;
  const REAL radiusA = 1;
  const REAL radiusB = 1;
  STRAND_MESH* strandMesh = new STRAND_MESH(restVertices, vertices,
                                                    E, nu, density, radiusA, radiusB);

  // create the integrator
  QUADRATIC_STRETCHING* stretchingEnergy = new QUADRATIC_STRETCHING(100.0);
  TIMESTEPPER* strandSolver = new TIMESTEPPER(*strandMesh, *stretchingEnergy);
  strandSolver->setDt(1.0 / 100.0);

  // constrain the first vertex
  VECTOR3 center0(-1.0,-0.2, 0.0);
  vector<KINEMATIC_SHAPE*> kinematicShapes;
  kinematicShapes.push_back(new CUBE(center0, 2.1));
  strandSolver->attachKinematicSurfaceConstraints(kinematicShapes[0]);

  // step forward nine timesteps
  for (int x = 0; x < 10; x++)
    strandSolver->solveQuasistatic(true);

  // get the current positions
  VECTOR positions = strandMesh->getPositions();
  VECTOR thetas = strandMesh->thetas();

  // initialize the known positions
  VECTOR positionsKnown(12);
  positionsKnown << -1, 0, 0, 0, 0, 0, 0.9994173123909491, 0.1243500915608951, -0.2291030282394728, 0.8139664087434525, 0.1410819371218592, -1.225833006238255;
  VECTOR thetasKnown(3);
  thetasKnown << 0, 0.7256850604191948, -1.252920931479943;

  // see if they match
  const REAL positionsDiff = (positionsKnown - positions).norm();
  const REAL thetasDiff = (thetasKnown - thetas).norm();

  cout << " position diff: " << positionsDiff << endl;
  cout << " theta diff: " << thetasDiff << endl; 
  if (positionsDiff + thetasDiff < 1e-8)
    cout << " Regression test PASSED" << endl;
  else
    cout << " Regression test FAILED" << endl;

  REQUIRE(positionsDiff < 1e-8);
  REQUIRE(thetasDiff < 1e-8);

  delete strandMesh;
  delete stretchingEnergy;
  delete strandSolver;
  delete kinematicShapes[0];
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void generateStrandGrid(vector<VECTOR3>& restVertices, vector<vector<int> >& strandIndices,
                        vector<KINEMATIC_SHAPE*>& kinematicShapes)
{
  const int _topStrands = 10;
  const REAL _topStrandSpacing = 0.2;
  const int _bottomStrands = 10;
  const REAL _bottomStrandSpacing = 0.2;

  // how far down the bottom strand is
  const REAL bottomHeight = 0.2;  // close up test
  const REAL topStart = -0.5 * _topStrandSpacing * _topStrands;
  const REAL bottomStart = 3.0 - 0.5 * _bottomStrandSpacing * _bottomStrands;

  const unsigned int totalTopPoints = 20;
  const unsigned int totalBottomPoints = 25;
  const REAL deltaTop = 3.0 / (REAL)(totalTopPoints - 1);
  const REAL deltaBottom = 3.0 / (REAL)(totalBottomPoints - 1);
  assert(totalTopPoints >= 4);

  for (unsigned int y = 0; y < _topStrands; y++)
  {
    for (unsigned int x = 0; x < totalTopPoints; x++)
      restVertices.push_back(VECTOR3(topStart  + y * _topStrandSpacing,0, deltaTop * x));

    // store the indices of individual strands
    vector<int> strand0;
    const unsigned int nextStart = y * totalTopPoints;
    //for (unsigned int x = 0; x < totalTopPoints; x++)
    for (unsigned int x = nextStart; x < nextStart + totalTopPoints; x++)
      strand0.push_back(x);
    strandIndices.push_back(strand0);
  } 

  for (unsigned int y = 0; y < _bottomStrands; y++)
  {
    const unsigned int nextStart = restVertices.size();
    for (unsigned int x = 0; x < totalBottomPoints; x++)
      restVertices.push_back(VECTOR3(-1.5 + x * deltaBottom,bottomHeight, bottomStart - y * _bottomStrandSpacing));
    
    // store the indices of individual strands
    vector<int> strand1;
    for (unsigned int x = nextStart ; x < nextStart + totalBottomPoints; x++)
      strand1.push_back(x);
    strandIndices.push_back(strand1);
  }
  
  const REAL cubeScale = 2.1;

  // constrain the first vertex
  for (unsigned int y = 0; y < _topStrands; y++)
  {
    VECTOR3 center1(topStart + y * _topStrandSpacing,0,3);
    kinematicShapes.push_back(new CUBE(center1, cubeScale * deltaTop));
  }
  for (unsigned int y = 0; y < _bottomStrands; y++)
  {
    VECTOR3 center2(-1.5,bottomHeight,bottomStart - y * _bottomStrandSpacing);
    kinematicShapes.push_back(new CUBE(center2, cubeScale * deltaBottom));
    VECTOR3 center3( 1.5,bottomHeight,bottomStart - y * _bottomStrandSpacing);
    kinematicShapes.push_back(new CUBE(center3, cubeScale * deltaBottom));
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
vector<int> runCollisionScene(STRAND_MESH& strandMesh, vector<KINEMATIC_SHAPE*>& kinematicShapes)
{
  // create the integrator
  STRAND::QUADRATIC_STRETCHING stretchingEnergy(strandMesh.kStretching());
  STRAND::TIMESTEPPER strandSolver(strandMesh, stretchingEnergy);
  strandSolver.setDt(1.0 / 300.0);

  strandSolver.hessianClampingEnabled() = true;
  strandSolver.collisionStiffness() = 10;
  strandSolver.collisionDampingBeta() = 0;

  for (unsigned int x = 0; x < kinematicShapes.size(); x++)
    strandSolver.attachKinematicSurfaceConstraints(kinematicShapes[x]);

  vector<int> totalCollisions;
  const VECTOR3 gravity = VECTOR3(0,981,0);
  //for (unsigned int x = 0; x < 100; x++)
  for (unsigned int x = 0; x < 20; x++)
  {
    strandSolver.externalForces().setZero();
    strandSolver.addGravity(gravity);
    strandSolver.solveDynamics(false);

    totalCollisions.push_back(strandMesh.totalCollisions());
  }

  return totalCollisions;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void fasterStrandCollisionDetection()
{
  cout << "==================================================================== " << endl;
  cout << " Matching STRAND_MESH_FASTER and STRAND_MESH collision detection" << endl;
  cout << "==================================================================== " << endl;
  // procedurally generate a scene
  std::vector<VECTOR3> restVertices;
  vector<vector<int> > strandIndices;
  vector<KINEMATIC_SHAPE*> kinematicShapes;
  generateStrandGrid(restVertices, strandIndices, kinematicShapes);

  const REAL E = 1e5; // things seem to misbehave here
  const REAL nu = 0.36;
  const REAL density = 1.32;
  const REAL radiusA = 0.005;
  const REAL radiusB = 0.005;
  const REAL collisionEps = 0.1;

  // get ground truth data
  STRAND_MESH strandMesh(restVertices, strandIndices,
                          E, nu, density, radiusA, radiusB);
  strandMesh.setCollisionEps(collisionEps);
  strandMesh.bendingForceFilterEnabled() = false;
  vector<int> groundTruthCollisions = runCollisionScene(strandMesh, kinematicShapes);

  // get accelerated version
  STRAND_MESH_FASTER strandMeshFaster(restVertices, strandIndices,
                                      E, nu, density, radiusA, radiusB);
  strandMeshFaster.setCollisionEps(collisionEps);
  strandMeshFaster.bendingForceFilterEnabled() = false;
  vector<int> fasterCollisions = runCollisionScene(strandMeshFaster, kinematicShapes);

  assert(groundTruthCollisions.size() == fasterCollisions.size());

  // count the differences
  int diff = 0;
  int totalCollisions;
  cout << " Collision differences: " << endl;
  for (unsigned int x = 0; x < fasterCollisions.size(); x++)
  {
    REAL currentDiff = abs(fasterCollisions[x] - groundTruthCollisions[x]);
    cout << currentDiff << " ";
    diff += currentDiff;

    totalCollisions += groundTruthCollisions[x];
  }
  cout << endl;
  cout << " Total collisions seen: " << totalCollisions << endl;

  REQUIRE(diff == 0);

  for (unsigned int x = 0; x < kinematicShapes.size(); x++)
    delete kinematicShapes[x];
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the force of an isotropic bending model
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestSpatialGradient(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending force for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 E0;
  E0.col(0) = vertices[0] - vertices[1];
  E0.col(1) = vertices[2] - vertices[1];

  REAL psi0 = material->psi(E0, theta0);
  MATRIX3x2 p = material->PK1(E0, theta0);

  MATRIX3 I = MATRIX3::Identity();
  MATRIX dEdx(6,9);
  dEdx.setZero();
  dEdx.block<3,3>(0,0) =  I;
  dEdx.block<3,3>(0,3) = -I;
  dEdx.block<3,3>(3,3) = -I;
  dEdx.block<3,3>(3,6) =  I;
  VECTOR9 gradientDirect = dEdx.transpose() * flatten(p);

  vector<VECTOR3> vertices0 = vertices;
  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR9 gradientNumerical;
    gradientNumerical.setZero();

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
        vector<VECTOR3> vertices1 = vertices0;
        vertices1[y][x] += eps;
  
        MATRIX3x2 E1;
        E1.col(0) = vertices1[0] - vertices1[1];
        E1.col(1) = vertices1[2] - vertices1[1];

        // get the new psi
        double psi = material->psi(E1, theta0);

        // store the finite difference
        gradientNumerical[3 * y + x] = (psi - psi0) / eps;
      }

    VECTOR9 diff = gradientDirect - gradientNumerical;
    REAL diffNorm = (fabs(diff.norm() / gradientDirect.norm())) / 9.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " direct: " << endl << gradientDirect << endl;
      cout << " finite diff: " << endl << gradientNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the force gradient of an isotropic bending model
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestSpatialHessian(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending force gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 E0;
  E0.col(0) = vertices[0] - vertices[1];
  E0.col(1) = vertices[2] - vertices[1];

  MATRIX3x2 p0 = material->PK1(E0, theta0);
  MATRIX6 H = material->hessian(E0, theta0);

  MATRIX3 I = MATRIX3::Identity();
  MATRIX dEdx(6,9);
  dEdx.setZero();
  dEdx.block<3,3>(0,0) =  I;
  dEdx.block<3,3>(0,3) = -I;
  dEdx.block<3,3>(3,3) = -I;
  dEdx.block<3,3>(3,6) =  I;
  MATRIX9 hessianDirect = dEdx.transpose() * H * dEdx;

  VECTOR9 f0 = dEdx.transpose() * flatten(p0);

  vector<VECTOR3> vertices0 = vertices;
  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX9 hessianNumerical;
    hessianNumerical.setZero();

    // for each of the degrees of the freedom
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
        vector<VECTOR3> vertices1 = vertices0;
        vertices1[y][x] += eps;
  
        MATRIX3x2 E1;
        E1.col(0) = vertices1[0] - vertices1[1];
        E1.col(1) = vertices1[2] - vertices1[1];

        MATRIX3x2 p1 = material->PK1(E1, theta0);
        VECTOR9 f1 = dEdx.transpose() * flatten(p1);

        // store the finite difference
        hessianNumerical.col(3 * y + x) = (f1 - f0) / eps;
      }

    MATRIX9 diff = hessianDirect - hessianNumerical;
    REAL diffNorm = (fabs(diff.norm() / hessianDirect.norm())) / 81.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " direct: " << endl << hessianDirect << endl;
      cout << " finite diff: " << endl << hessianNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the force gradient of an isotropic bending model
//////////////////////////////////////////////////////////////////////////////
bool centeredTestSpatialHessian(const HOBAK::STRAND::ISOTROPIC_BENDING* material, const vector<VECTOR3>& vertices)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Isotropic Bending centered force gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 E0;
  E0.col(0) = vertices[0] - vertices[1];
  E0.col(1) = vertices[2] - vertices[1];

  //MATRIX3x2 p0 = material->PK1(E0, theta0);
  //REAL p0 = material->psi(E0, theta0);
  MATRIX6 H = material->hessian(E0, theta0);

  MATRIX3 I = MATRIX3::Identity();
  MATRIX dEdx(6,9);
  dEdx.setZero();
  dEdx.block<3,3>(0,0) =  I;
  dEdx.block<3,3>(0,3) = -I;
  dEdx.block<3,3>(3,3) = -I;
  dEdx.block<3,3>(3,6) =  I;
  MATRIX9 hessianDirect = dEdx.transpose() * H * dEdx;

  const vector<VECTOR3> vertices0 = vertices;
  double eps = 1e-2;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-5)
  {
    MATRIX9 hessianNumerical;
    hessianNumerical.setZero();

    // for each of the degrees of the freedom
    for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
        vector<VECTOR3> verticesPlus = vertices0;
        verticesPlus[j][i] += eps;
        MATRIX3x2 EPlus;
        EPlus.col(0) = verticesPlus[0] - verticesPlus[1];
        EPlus.col(1) = verticesPlus[2] - verticesPlus[1];
        
        vector<VECTOR3> verticesMinus = vertices0;
        verticesMinus[j][i] -= eps;
        MATRIX3x2 EMinus;
        EMinus.col(0) = verticesMinus[0] - verticesMinus[1];
        EMinus.col(1) = verticesMinus[2] - verticesMinus[1];

        for (int y = 0; y < 3; y++)
          for (int x = 0; x < 3; x++)
          {
            vector<VECTOR3> verticesPP = verticesPlus;
            verticesPP[y][x] += eps;

            MATRIX3x2 EPP;
            EPP.col(0) = verticesPP[0] - verticesPP[1];
            EPP.col(1) = verticesPP[2] - verticesPP[1];
            const REAL pp = material->psi(EPP, theta0);

            vector<VECTOR3> verticesPM = verticesPlus;
            verticesPM[y][x] -= eps;

            MATRIX3x2 EPM;
            EPM.col(0) = verticesPM[0] - verticesPM[1];
            EPM.col(1) = verticesPM[2] - verticesPM[1];
            const REAL pm = material->psi(EPM, theta0);

            vector<VECTOR3> verticesMP = verticesMinus;
            verticesMP[y][x] += eps;

            MATRIX3x2 EMP;
            EMP.col(0) = verticesMP[0] - verticesMP[1];
            EMP.col(1) = verticesMP[2] - verticesMP[1];
            const REAL mp = material->psi(EMP, theta0);

            vector<VECTOR3> verticesMM = verticesMinus;
            verticesMM[y][x] -= eps;

            MATRIX3x2 EMM;
            EMM.col(0) = verticesMM[0] - verticesMM[1];
            EMM.col(1) = verticesMM[2] - verticesMM[1];
            const REAL mm = material->psi(EMM, theta0);

            /*
            if (x == 0 && y == 0 && i == 0 && j == 0)
            {
              cout << " pp: " << pp << endl;
              cout << " mp: " << mp << endl;
              cout << " pm: " << pm << endl;
              cout << " mm: " << mm << endl;
              cout << " finite diff: " << (pp + mm - mp - pm) / (4.0 * eps * eps) << endl;
            }
            */

            // use the star-shaped mixed finite-difference
            hessianNumerical(3 * j + i, 3 * y + x) = (pp + mm - mp - pm) / (4.0 * eps * eps);
          }
      }

    MATRIX9 diff = hessianDirect - hessianNumerical;
    REAL diffNorm = (fabs(diff.norm() / hessianDirect.norm())) / 81.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (minSeen > 1e-5)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " direct: " << endl << hessianDirect << endl;
      cout << " finite diff: " << endl << hessianNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      return false;
    }
    else
      eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-5)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the force gradient of an isotropic bending model
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTwistingSpatialHessian(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Twisting centered force gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  vector<VECTOR3> vertices0(4);
  for (int x = 0; x < 4; x++)
    vertices0[x] = randomVector3();

#if 0
  MATRIX3 F0 = computeFtwist(vertices0);
  // // rank 4 correction
  // MATRIX9x12 pFpX = computeTwistPFPxOrig(vertices0);

  // MATRIX9 H = material->hessian(F0);
  // MATRIX12 hessianDirect = pFpX.transpose() * H * pFpX + computeExtraTerm(vertices0, material);

  // // testing analytical eigen decomposition of rank 4 correction
  // MATRIX12 hessExtra = computeExtraTerm(vertices0, material);
  // MATRIX12 hessExtraClamped = computeExtraTermClamped(vertices0, material);
  // if((hessExtra - hessExtraClamped).norm()/hessExtra.norm()>1e-5)
  // {
  //   cout << " TEST FAILED!!!!!" << endl;
  //   cout << " Original: " << endl << hessExtra << endl;
  //   cout << " clamped: " << endl << hessExtraClamped << endl;
  // }
  // MATRIX12 hessianDiffDirect = computeExtraTerm(vertices0, material);
  MATRIX9x12 pFpX = computeTwistPFPx(vertices0);
  MATRIX9x12 delpFpX = computeDelTwistPFPx(vertices0);

  MATRIX9 H = material->hessian(F0);
  MATRIX12 hessianDirect = pFpX.transpose() * H * pFpX - delpFpX.transpose() * H * delpFpX;
#else
  MATRIX12 hessianDirect = material->forceGradient(vertices0);
#endif

  double eps = 1e-2;
  int e = 0;
  double minSeen = FLT_MAX;
  const int maxTries = 3;
  while (e <= maxTries)
  {
    MATRIX12 hessianNumerical;
    hessianNumerical.setZero();

    // for each of the degrees of the freedom
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++)
      {
        vector<VECTOR3> verticesPlus = vertices0;
        verticesPlus[j][i] += eps;
        
        vector<VECTOR3> verticesMinus = vertices0;
        verticesMinus[j][i] -= eps;

        for (int y = 0; y < 4; y++)
          for (int x = 0; x < 3; x++)
          {
            vector<VECTOR3> verticesPP = verticesPlus;
            verticesPP[y][x] += eps;
            const MATRIX3 FPP = computeFtwist(verticesPP);
            const REAL pp = material->psi(FPP);

            vector<VECTOR3> verticesPM = verticesPlus;
            verticesPM[y][x] -= eps;
            const MATRIX3 FPM = computeFtwist(verticesPM);
            const REAL pm = material->psi(FPM);

            vector<VECTOR3> verticesMP = verticesMinus;
            verticesMP[y][x] += eps;
            const MATRIX3 FMP = computeFtwist(verticesMP);
            const REAL mp = material->psi(FMP);

            vector<VECTOR3> verticesMM = verticesMinus;
            verticesMM[y][x] -= eps;
            const MATRIX3 FMM = computeFtwist(verticesMM);
            const REAL mm = material->psi(FMM);

            // use the star-shaped mixed finite-difference
            hessianNumerical(3 * j + i, 3 * y + x) = (pp + mm - mp - pm) / (4.0 * eps * eps);
          }
      }

    MATRIX12 diff = hessianDirect - hessianNumerical;
    // MATRIX12 diffdiff = diff + hessianDiffDirect;
    REAL diffNorm = (fabs(diff.norm() / hessianDirect.norm())) / 144.0;
    // REAL diffNorm = (fabs(diffdiff.norm() / diff.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (minSeen > 1e-6 && e == maxTries)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " direct: " << endl << hessianDirect << endl;
      cout << " finite diff: " << endl << hessianNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      // cout << " direct: " << endl << hessianDiffDirect << endl;
      // cout << " finite diff: " << endl << diff << endl;
      // cout << " error: " << endl << diffdiff << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on the force gradient of an isotropic bending model
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestTwistingSpatialHessianRank4(const HOBAK::VOLUME::TET_STRAND_TWIST* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Strand Twisting centered force gradient for " << material->name().c_str() <<  " with rank-4 correction" << endl;
  cout << "=============================================================== " << endl;

  vector<VECTOR3> vertices0(4);
  for (int x = 0; x < 4; x++)
    vertices0[x] = randomVector3();

#if 0
  MATRIX3 F0 = computeFtwist(vertices0);
  // // rank 4 correction
  MATRIX9x12 pFpX = computeTwistPFPxOrig(vertices0);

  MATRIX9 H = material->hessian(F0);
  MATRIX12 hessianDirect = pFpX.transpose() * H * pFpX + computeExtraTerm(vertices0, material);
#else
  MATRIX12 hessianDirect = material->forceGradientRank4(vertices0);
#endif

  double eps = 1e-2;
  int e = 0;
  double minSeen = FLT_MAX;
  const int maxTries = 3;
  while (e <= maxTries)
  {
    MATRIX12 hessianNumerical;
    hessianNumerical.setZero();

    // for each of the degrees of the freedom
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++)
      {
        vector<VECTOR3> verticesPlus = vertices0;
        verticesPlus[j][i] += eps;
        
        vector<VECTOR3> verticesMinus = vertices0;
        verticesMinus[j][i] -= eps;

        for (int y = 0; y < 4; y++)
          for (int x = 0; x < 3; x++)
          {
            vector<VECTOR3> verticesPP = verticesPlus;
            verticesPP[y][x] += eps;
            const MATRIX3 FPP = computeFtwist(verticesPP);
            const REAL pp = material->psi(FPP);

            vector<VECTOR3> verticesPM = verticesPlus;
            verticesPM[y][x] -= eps;
            const MATRIX3 FPM = computeFtwist(verticesPM);
            const REAL pm = material->psi(FPM);

            vector<VECTOR3> verticesMP = verticesMinus;
            verticesMP[y][x] += eps;
            const MATRIX3 FMP = computeFtwist(verticesMP);
            const REAL mp = material->psi(FMP);

            vector<VECTOR3> verticesMM = verticesMinus;
            verticesMM[y][x] -= eps;
            const MATRIX3 FMM = computeFtwist(verticesMM);
            const REAL mm = material->psi(FMM);

            // use the star-shaped mixed finite-difference
            hessianNumerical(3 * j + i, 3 * y + x) = (pp + mm - mp - pm) / (4.0 * eps * eps);
          }
      }

    MATRIX12 diff = hessianDirect - hessianNumerical;
    // MATRIX12 diffdiff = diff + hessianDiffDirect;
    REAL diffNorm = (fabs(diff.norm() / hessianDirect.norm())) / 144.0;
    // REAL diffNorm = (fabs(diffdiff.norm() / diff.norm())) / 144.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (minSeen > 1e-6 && e == maxTries)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " direct: " << endl << hessianDirect << endl;
      cout << " finite diff: " << endl << hessianNumerical << endl;
      cout << " diff: " << endl << diff << endl;
      // cout << " direct: " << endl << hessianDiffDirect << endl;
      // cout << " finite diff: " << endl << diff << endl;
      // cout << " error: " << endl << diffdiff << endl;
      return false;
    }
    eps *= 0.1;
    e++;
  }
  if (minSeen < 1e-6)
    cout << " TEST PASSED. " << endl;
  else
  {
    cout << " TEST FAILED. " << endl;
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// test out strand energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Strand Energy tests", "[strand tests]" ) 
{
  using namespace HOBAK::STRAND;
  VECTOR3 f = randomVector3(2.0);

  vector<VECTOR3> p(2);
  p[0] = randomVector3(2.0);
  p[1] = randomVector3(2.0);

  SECTION("Testing Quadratic Stretching") 
  {
    QUADRATIC_STRETCHING quadratic(123.0);
    REQUIRE(convergenceTestPK1(&quadratic, f) == true);
    REQUIRE(convergenceTestHessian(&quadratic, f) == true);
    REQUIRE(testClampedHessian(&quadratic) == true);
    REQUIRE(convergenceTestGradient(&quadratic, p) == true);
    REQUIRE(convergenceTestSpatialHessian(&quadratic, p) == true);
    REQUIRE(testSpatialClampedHessian(&quadratic) == true);
  }
  
  SECTION("Testing Theta Bending") 
  {
    MATRIX3x2 E = randomMatrix3x2();
    ISOTROPIC_THETA theta(123.0);
    REQUIRE(convergenceTestPK1(&theta, E) == true);
    REQUIRE(convergenceTestHessian(&theta, E) == true);
    REQUIRE(testClampedHessian(&theta) == true);
  }

  /*
  {
    MATRIX3x2 E;
    E << 1,1,1,0,0,0;
    cout << "E: " << endl << E << endl; 
    QUADRATIC_BENDING quadratic(1.0);
    QUADRATIC_UNIT_BENDING unit(1.0);

    cout << " quadratic: " << endl << quadratic.hessian(E,M_PI * 0.5) << endl;
    cout << " unit: " << endl << unit.hessian(E,M_PI * 0.5) << endl;

  }
  */

  SECTION("Testing Quadratic Bending") 
  {
    MATRIX3x2 E = randomMatrix3x2();
    QUADRATIC_BENDING quadratic(123.0);
    REQUIRE(convergenceTestPK1(&quadratic, E) == true);
    REQUIRE(convergenceTestHessian(&quadratic, E) == true);
    
    REQUIRE(testPathological() == true);
    REQUIRE(testClampedHessian(&quadratic) == true);
  }

  /*
  SECTION("Small strand regression test")
  {
    strandRegressionScene();
  } 
  */
  SECTION("Verify collision detection in STRAND_MESH_FASTER matched STRAND_MESH")
  {
    fasterStrandCollisionDetection();
  }
  SECTION("Verify the STRAND_MESH_FASTER didn't break anything")
  {
    fasterStrandMeshTest();
  }

  SECTION("Convergence test sine bending")
  {
    const STRAND::SIN_BENDING bending;
    REQUIRE(convergenceTestEdgeGradientSinBending(bending) == true);
    REQUIRE(convergenceTestThetaGradientSinBending(bending) == true);
    REQUIRE(convergenceTestGradientSinBending(bending) == true);
    REQUIRE(convergenceTestEdgeHessianSinBending(bending) == true);
    REQUIRE(convergenceTestHessianSinBending(bending) == true);
  }
  SECTION("Convergence test half bending")
  {
    const STRAND::HALF_BENDING bending; 
    REQUIRE(convergenceTestEdgeGradientSinBending(bending) == true);
    REQUIRE(convergenceTestEdgeHessianSinBending(bending) == true);
    REQUIRE(convergenceTestThetaGradientSinBending(bending) == true);
    REQUIRE(convergenceTestGradientSinBending(bending) == true);
    REQUIRE(convergenceTestHessianSinBending(bending) == true);
  }
  /*
  SECTION("Testing Strand twisting against a known answer")
  {
    testTwisting();
  } 
  SECTION("Testing DVT energy")
  {
    convergenceTestGradientDVT();
    convergenceTestHessianDVT();
  }
  */

  SECTION("Testing Quadratic Unit Bending")
  {
    MATRIX3x2 E = randomMatrix3x2();
    QUADRATIC_UNIT_BENDING quadratic(123.0);
    REQUIRE(convergenceTestPK1(&quadratic, E) == true);
    REQUIRE(convergenceTestHessian(&quadratic, E) == true);
    REQUIRE(convergenceTestPK1Inverted(&quadratic, E) == true);
    REQUIRE(convergenceTestHessianInverted(&quadratic, E) == true);
    
    vector<VECTOR3> vertices(3);
    for (int x = 0; x < 3; x++)
      vertices[x] = randomVector3();

    REQUIRE(convergenceTestSpatialGradient(&quadratic, vertices) == true);
    REQUIRE(convergenceTestSpatialHessian(&quadratic, vertices) == true);
    REQUIRE(centeredTestSpatialHessian(&quadratic, vertices) == true);
  }

  SECTION("Testing TAN Bending")
  {
    MATRIX3x2 E = randomMatrix3x2();
    TAN_BENDING tanBend(123.0);
    REQUIRE(convergenceTestPK1(&tanBend, E) == true);
    REQUIRE(convergenceTestHessian(&tanBend, E) == true);
    REQUIRE(convergenceTestPK1Inverted(&tanBend, E) == true);
    REQUIRE(convergenceTestHessianInverted(&tanBend, E) == true);
    
    vector<VECTOR3> vertices(3);
    for (int x = 0; x < 3; x++)
      vertices[x] = randomVector3();

    REQUIRE(convergenceTestSpatialGradient(&tanBend, vertices) == true);
    REQUIRE(convergenceTestSpatialHessian(&tanBend, vertices) == true);
    // REQUIRE(centeredTestSpatialHessian(&tanBend, vertices) == true);
  }

  SECTION("Testing Quadratic ANOTHER Bending")
  {
    MATRIX3x2 E = randomMatrix3x2();
    QUADRATIC_F_BENDING quadBend(123.0);
    REQUIRE(convergenceTestPK1(&quadBend, E) == true);
    REQUIRE(convergenceTestHessian(&quadBend, E) == true);
    REQUIRE(convergenceTestPK1Inverted(&quadBend, E) == true);
    REQUIRE(convergenceTestHessianInverted(&quadBend, E) == true);
    
    vector<VECTOR3> vertices(3);
    for (int x = 0; x < 3; x++)
      vertices[x] = randomVector3();

    REQUIRE(convergenceTestSpatialGradient(&quadBend, vertices) == true);
    REQUIRE(convergenceTestSpatialHessian(&quadBend, vertices) == true);
    // REQUIRE(centeredTestSpatialHessian(&quadBend, vertices) == true);
  }


  // SECTION("Convergence test Quadratic Angle Spring Twisting")
  // {
  //   MATRIX3 F;
  //   VECTOR3 e0 = randomVector3();
  //   VECTOR3 e1 = randomVector3();
  //   VECTOR3 e2 = randomVector3();
  //   VECTOR3 e1Hat = e1.normalized();
  //   e0 = e0 - e0.dot(e1Hat) * e1Hat;
  //   e2 = e2 - e2.dot(e1Hat) * e1Hat;
  //   F.col(0) = e0;
  //   F.col(1) = e2;
  //   F.col(2) = e1;

  //   HOBAK::VOLUME::TET_STRAND_TWIST twisting(10.0,M_PI/3.0);
  //   REQUIRE(convergenceTestPK1(&twisting, F));
  //   REQUIRE(convergenceTestTwistHessian(&twisting, F));
  //   REQUIRE(convergenceTestTwistForcesConst(twisting));
  //   REQUIRE(convergenceTestTwistingSpatialGradient(&twisting));
  //   REQUIRE(convergenceTestTwistingSpatialHessian(&twisting));
  //   REQUIRE(convergenceTestTwistingSpatialHessianRank4(&twisting));
  //   REQUIRE(testClampedForceGradient(&twisting));
  //   //REQUIRE(testClampedForceGradientRank4(&twisting));
  // }
}
