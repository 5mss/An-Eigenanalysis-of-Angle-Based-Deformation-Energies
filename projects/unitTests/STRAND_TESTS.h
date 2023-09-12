using namespace HOBAK;
using namespace std;

const REAL theta0 = M_PI * 0.5;



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
// test out strand energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Strand Energy tests", "[strand tests]" ) 
{
  using namespace HOBAK::STRAND;
  VECTOR3 f = randomVector3(2.0);

  vector<VECTOR3> p(2);
  p[0] = randomVector3(2.0);
  p[1] = randomVector3(2.0);



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

  SECTION("Testing Quadratic F Bending")
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
