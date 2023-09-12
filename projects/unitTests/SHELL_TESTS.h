using namespace HOBAK;
using namespace std;


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
// test out isotropic energies
//////////////////////////////////////////////////////////////////////////////
TEST_CASE( "Isotropic Shell Energy tests", "[isotropic shells]" ) 
{
  using namespace HOBAK::SHELL;
  MATRIX3x2 F = randomMatrix3x2();

  
  vector<VECTOR3> flap;
  flap.push_back(VECTOR3(-1,-1,0));
  flap.push_back(VECTOR3(0,0,-1));
  flap.push_back(VECTOR3(0,0,1));
  flap.push_back(VECTOR3(1,-1,0));


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
  }
}
