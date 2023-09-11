using namespace HOBAK;
using namespace std;

//////////////////////////////////////////////////////////////////////////////
// do a convergence test on Hessian
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const HOBAK::SHELL::STRETCHING* material, const MATRIX3x2 &F)
{
  using namespace HOBAK;
  using namespace std;

  cout << "=============================================================== " << endl;
  cout << " VERIFYING Shell Hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  MATRIX3x2 U;
  MATRIX2 V;
  VECTOR2 Sigma;
  svd(F, U, Sigma, V);
  MATRIX6 dPdF = material->hessian(F);
  MATRIX3x2 P = material->PK1(F);

  double eps = 1e-4;
  int e = 0;
  double minSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX6 finiteDiff;
    int column = 0;

    // for each of the degrees of the freedom
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 3; x++, column++)
      {
        MATRIX3x2 Fnew = F;
        Fnew(x,y) += eps;

        // get the new psi
        MATRIX3x2 Pnew = material->PK1(Fnew);

        // store the finite difference
        MATRIX3x2 diff = (Pnew - P) / eps;
        finiteDiff.col(column) = flatten(diff);
      }

    MATRIX6 diff = dPdF - finiteDiff;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 12.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    MATRIX6 div = finiteDiff;
    for (int y = 0; y < 6; y++)
      for (int x = 0; x < 6; x++)
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
// do a convergence test on the bending energy hessian
/////////////////////////////////////////////////////////////////////////////
bool convergenceTestHessian(const HOBAK::SHELL::BENDING* material, 
                            const vector<VECTOR3>& flap)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Shell Bending hessian for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  //const REAL theta = material->restAngle(flap);
  const REAL theta = 1.0;
  // cout<<"--------------Psi: "<<material->psi(flap, theta)<<endl;
  VECTOR12 g0 = material->gradient(flap, theta);
  MATRIX12 H = material->hessian(flap, theta);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    MATRIX12 finiteDiff;

    // for each of the degrees of the freedom
    int index = 0;
    for (int y = 0; y < 4; y++)
      for (int x = 0; x < 3; x++, index++)
      {
        vector<VECTOR3> flapNew = flap;
        flapNew[y][x] += eps;

        // get the new psi
        VECTOR12 g = material->gradient(flapNew, 1.0);

        // store the finite difference
        finiteDiff.col(index) = (g - g0) / eps;
      }

    MATRIX12 diff = H - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / H.norm())) / 144.0;
    REAL absDiffNorm = fabs(diff.norm()) / 144.0;
    if (relativeDiffNorm < relativeMinSeen)
      relativeMinSeen = relativeDiffNorm;
    if (absDiffNorm < absMinSeen)
      absMinSeen = absDiffNorm;
    cout << "eps: " << eps << " relative diff: " << relativeDiffNorm << "\t abs diff: " << absDiffNorm << endl;

    if (e == 4 && relativeMinSeen > 1e-6 && absDiffNorm > 1e-8)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " H: " << endl << clampSmalls(H) << endl;
      cout << " finite diff: " << endl << clampSmalls(finiteDiff) << endl;
      cout << " diff: " << endl << clampSmalls(diff) << endl;
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
// do a convergence test on the bending energy gradient
//////////////////////////////////////////////////////////////////////////////
bool convergenceTestGradient(const HOBAK::SHELL::BENDING* material, 
                             const vector<VECTOR3>& flap)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Shell Bending gradient for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  const REAL theta = material->restAngle(flap);
  // const REAL theta = 1.0;

  REAL psi0 = material->psi(flap, theta);
  VECTOR12 g = material->gradient(flap, theta);

  double eps = 1e-4;
  int e = 0;
  double relativeMinSeen = FLT_MAX;
  double absMinSeen = FLT_MAX;
  while (eps > 1e-8)
  {
    VECTOR12 finiteDiff;

    // for each of the degrees of the freedom
    int index = 0;
    for (int y = 0; y < 4; y++)
      for (int x = 0; x < 3; x++, index++)
      {
        vector<VECTOR3> flapNew = flap;
        flapNew[y][x] += eps;

        // get the new psi
        double psi = material->psi(flapNew, theta);

        // store the finite difference
        finiteDiff[index] = (psi - psi0) / eps;

        /*
        if (index == 0)
        {
          cout << " psi0: " << psi0 << "\t psi: " << psi << endl;
          cout << " finite diff: " << (psi - psi0) / eps << endl;
          cout << " diff: " << psi - psi0 << endl;
          cout << " diff 2 pi: " << (psi - psi0) + (2.0 * M_PI) << endl;
        }
        */
      }

    VECTOR12 diff = g - finiteDiff;
    REAL relativeDiffNorm = (fabs(diff.norm() / g.norm())) / 12.0;

    // need this check in case the vector is just all zeros
    REAL absDiffNorm = fabs(diff.norm()) / 12.0;
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
bool convergenceTestPK1(const HOBAK::SHELL::STRETCHING* material, const MATRIX3x2 &F)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Shell P for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;

  REAL psi0 = material->psi(F);
  MATRIX3x2 P = material->PK1(F);

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
        MATRIX3x2 Fnew = F;
        Fnew(x,y) += eps;

        // get the new psi
        double psi = material->psi(Fnew);

        // store the finite difference
        finiteDiffP(x, y) = (psi - psi0) / eps;
      }

    MATRIX3x2 diff = P - finiteDiffP;
    REAL diffNorm = (fabs(diff.norm() / P.norm())) / 6.0;
    if (diffNorm < minSeen)
      minSeen = diffNorm;
    cout << "eps: " << eps << " diff: " << diffNorm << endl;

    if (e == 4 && minSeen > 1e-6)
    {
      cout << " TEST FAILED!!!!!" << endl;
      cout << " P: " << endl << P << endl;
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
bool testClampedHessian(const HOBAK::SHELL::STRETCHING* material)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  MATRIX3x2 F = randomMatrix3x2();

  MATRIX6 original, clamped, diff;
  REAL diffNorm;

  // verify that with all positive, we still get the same Hessian
  original = clampEigenvalues(material->hessian(F));
  clamped  = material->clampedHessian(F);
  diff = original - clamped;
  diffNorm = fabs(diff.norm() / original.norm());
  cout << " Diff: " << diffNorm << "\t";

  if (diffNorm < 1e-8 || original.norm() < 1e-8)
    cout << " All positive Sigma test: PASSED " << endl;
  else
  {
    cout << " All positive Sigma test: FAILED " << endl;
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
  
  // if it's not invertible, don't pass it an inverted F
  const REAL leftBound = -1.0;

  // try lots of Fs here, as this case tends to be tricky
  std::mt19937 gen(314159);
  std::uniform_real_distribution<REAL> dist(leftBound, 1.0);
  for (int x = 0; x < 20; x++)
  {
    // try the all-under-compression case, see if it still holds 
    //VECTOR2 sigmas = VECTOR2(dist(gen), dist(gen));
    MATRIX3x2 sigmas;
    sigmas.setZero();
    sigmas(0,0) = dist(gen);
    sigmas(1,1) = dist(gen);
    original = clampEigenvalues(material->hessian(sigmas));
    clamped  = material->clampedHessian(sigmas);
    diff = original - clamped;
    diffNorm = fabs(diff.norm() / original.norm());
    cout << " Diff: " << diffNorm << "\t";

    if (diffNorm < 1e-8 || original.norm() < 1e-8)
      cout << " Compression test: PASSED " << endl;
    else
    {
      MATRIX6 unclamped = material->hessian(sigmas);
      cout << " Compression test: FAILED " << endl;
      cout << " Unclamped: " << endl << unclamped << endl;
      cout << " Numerical clamped: " << endl << original << endl;
      cout << " Analytic clamped: " << endl << clamped << endl;
      cout << " Diff: " << endl << diff << endl;

      cout << " Numerical eigs: " << eigenvalues(original).transpose() << endl;
      cout << " Analytic eigs:  " << eigenvalues(clamped).transpose() << endl;
      cout << " Unclamped eigs: " << eigenvalues(unclamped).transpose() << endl;
      cout << " Sigmas:         " << sigmas.transpose() << endl;
      return false;
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// see if eigenvalue clamping works
//////////////////////////////////////////////////////////////////////////////
bool testClampedHessian(const HOBAK::SHELL::BENDING* material, const vector<VECTOR3>& flap)
{
  cout << "=============================================================== " << endl;
  cout << " VERIFYING Eigenvalue clamping for " << material->name().c_str() << endl;
  cout << "=============================================================== " << endl;
  const REAL theta = 1.0;
  // cout<<"--------------Psi: "<<material->psi(flap, theta)<<endl;
  MATRIX12 H = material->hessian(flap, theta);
  MATRIX12 Heig = material->clampedHessian(flap, theta);
  if((H - Heig).norm()/H.norm() > 1e-6){
    cout << " Eigen decomposition test: FAILED " << endl;
    cout << " H: " << endl << H << endl;
    cout << " Heig: " << endl << Heig << endl;
    return false;
  }
  cout<<"Bending Eigen decomposition PASSED.\n";
  return true;
}

//////////////////////////////////////////////////////////////////////////////
// test out isotropic energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Isotropic Shell Energy tests", "[isotropic shells]" ) 
{
  using namespace HOBAK::SHELL;
  MATRIX3x2 F = randomMatrix3x2();

  SECTION("Testing ARAP") 
  {
    ARAP arap(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&arap, F) == true);
    REQUIRE(convergenceTestHessian(&arap, F) == true);
    REQUIRE(testClampedHessian(&arap) == true);
  }
  SECTION("Testing StVK")
  {
    STVK stvk(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&stvk, F) == true);
    REQUIRE(convergenceTestHessian(&stvk, F) == true);
    REQUIRE(testClampedHessian(&stvk) == true);
  }
  SECTION("Testing Baraff-Witkin Stretch")
  {
    BW_STRETCH bw(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&bw, F) == true);
    REQUIRE(convergenceTestHessian(&bw, F) == true);
    REQUIRE(testClampedHessian(&bw) == true);
  }

  SECTION("Testing Baraff-Witkin Shear")
  {
    BW_SHEAR bw(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&bw, F) == true);
    REQUIRE(convergenceTestHessian(&bw, F) == true);
    REQUIRE(testClampedHessian(&bw) == true);
  }
  
  SECTION("Testing Baraff-Witkin")
  {
    BARAFF_WITKIN bw(1.0, 1.0);
    REQUIRE(convergenceTestPK1(&bw, F) == true);
    REQUIRE(convergenceTestHessian(&bw, F) == true);
    //REQUIRE(testClampedHessian(&bw) == true);
  }

  vector<VECTOR3> flap;
  flap.push_back(VECTOR3(-1,-1,0));
  flap.push_back(VECTOR3(0,0,-1));
  flap.push_back(VECTOR3(0,0,1));
  flap.push_back(VECTOR3(1,-1,0));

  SECTION("Testing Bending angle computation")
  {
    cout << "========================================================================== " << endl;
    cout << " TESTING Shell bending angle computation" << endl;
    cout << "========================================================================== " << endl;
    vector<VECTOR3> thetaFlap = flap;
    const VECTOR3 original = thetaFlap[3];

    REAL angle;
    THETA bending(1.0);
    
    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][0] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 135"<< endl;
    REQUIRE(fabs(angle - 135.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3] = thetaFlap[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 180"<< endl;
    //REQUIRE(fabs(angle - 180.0) < 1e-8);
    //REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    //REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3][0] = original[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 90"<< endl;
    REQUIRE(fabs(angle - 90.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3][1] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 45"<< endl;
    REQUIRE(fabs(angle - 45.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3] = -1.0 * thetaFlap[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 0"<< endl;
    REQUIRE(fabs(angle - 0.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3][0] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -45"<< endl;
    REQUIRE(fabs(angle + 45.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][1] *= -1;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -90"<< endl;
    REQUIRE(fabs(angle + 90.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][1] = 0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -135"<< endl;
    REQUIRE(fabs(angle + 135.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
  }

  SECTION("Testing faster Bending angle computation")
  {
    cout << "========================================================================== " << endl;
    cout << " TESTING faster Shell bending angle computation" << endl;
    cout << "========================================================================== " << endl;
    vector<VECTOR3> thetaFlap = flap;
    const VECTOR3 original = thetaFlap[3];

    REAL angle;
    THETA_FASTER bending(1.0);
    
    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][0] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 135"<< endl;
    REQUIRE(fabs(angle - 135.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3] = thetaFlap[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 180"<< endl;
    //REQUIRE(fabs(angle - 180.0) < 1e-8);
    //REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    //REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3][0] = original[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 90"<< endl;
    REQUIRE(fabs(angle - 90.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3][1] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 45"<< endl;
    REQUIRE(fabs(angle - 45.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3] = -1.0 * thetaFlap[0];
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be 0"<< endl;
    REQUIRE(fabs(angle - 0.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3][0] = 0.0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -45"<< endl;
    REQUIRE(fabs(angle + 45.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
    
    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][1] *= -1;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -90"<< endl;
    REQUIRE(fabs(angle + 90.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);

    thetaFlap[3] = thetaFlap[0];
    thetaFlap[3][1] = 0;
    angle = bending.restAngle(thetaFlap);
    angle *= 360.0 / (2.0 * M_PI);
    cout << "========================================================================== " << endl;
    cout << "========================================================================== " << endl;
    cout << "Angle " << angle << " should be -135"<< endl;
    REQUIRE(fabs(angle + 135.0) < 1e-8);
    REQUIRE(convergenceTestGradient(&bending, thetaFlap) == true);
    REQUIRE(convergenceTestHessian(&bending, thetaFlap) == true);
  }

  // SECTION("Testing Bending Spring")
  // {
  //   BENDING_SPRING bending(1.0);
  //   REQUIRE(convergenceTestGradient(&bending, flap) == true);
  //   REQUIRE(convergenceTestHessian(&bending, flap) == true);
  // }

  SECTION("Testing Dihedral Bending")
  {
    DIHEDRAL bending(1.0);
    REQUIRE(convergenceTestGradient(&bending, flap) == true);
    REQUIRE(convergenceTestHessian(&bending, flap) == true);
  }

  SECTION("Testing QUADRATIC_F_BENDING Bending")
  {
    QUADRATIC_F_BENDING bending(1.0);
    vector<VECTOR3> flapRand;
    flapRand.push_back(randomVector3());
    flapRand.push_back(randomVector3());
    flapRand.push_back(randomVector3());
    flapRand.push_back(randomVector3());
    REQUIRE(convergenceTestGradient(&bending, flap) == true);
    REQUIRE(convergenceTestHessian(&bending, flapRand) == true);
    flapRand[1] = VECTOR3(0.6488, 0.9653, 0.4605);
    flapRand[2] = VECTOR3(-0.3122, 0.1681, -0.7845);
    flapRand[0] = VECTOR3(0.8126, 0.7593, 0.6355);
    flapRand[3] = VECTOR3(-0.4785, 0.1887, -0.9550);
    REQUIRE(testClampedHessian(&bending, flapRand) == true);
  }
}
