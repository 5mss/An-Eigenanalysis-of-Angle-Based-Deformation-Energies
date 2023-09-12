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
#ifndef FIELD_2D_H
#define FIELD_2D_H

#include <cmath>
#include <string>
#include "SETTINGS.h"

// macro to debug 2D fields, assumes that fieldViewer has been built
// (not actually part of ANGLE yet)
#ifndef VARNAME
#define VARNAME(x) #x
#endif
#ifndef FIELDVIEW2D
#define FIELDVIEW2D(x) FIELD_2D::fieldViewer(x, VARNAME(x)); sleep(1);
#endif

namespace HOBAK {

class FIELD_2D {
public:
  FIELD_2D();
  FIELD_2D(const int& rows, const int& cols);
  FIELD_2D(const FIELD_2D& m);
  FIELD_2D(const MATRIX& m);
  ~FIELD_2D();

  // accessors
  inline REAL& operator()(int x, int y) { return _data[y * _xRes + x]; };
  const REAL operator()(int x, int y) const { return _data[y * _xRes + x]; };
  inline REAL& operator[](int x) { return _data[x]; };
  const REAL operator[](int x) const { return _data[x]; };
  REAL* data() { return _data; };
  const int xRes() const { return _xRes; };
  const int yRes() const { return _yRes; };
  const int totalCells() const { return _totalCells; };

  // a safe, slow, toroidal data accessor -- does all bounds checking for you
  inline REAL safe(int x, int y);

  // common field operations
  void clear();
  void normalize();
  FIELD_2D& abs();

  REAL min();
  REAL max();

  // field maximum cell index
  VECTOR3 maxIndex();

  // field minimum cell index
  VECTOR3 minIndex();

  // take the log
  void log(REAL base = 2.0);
 
  // IO functions
  void writePPM(std::string filename);
  void writeMatlab(std::string filename, std::string variableName) const;
  void write(std::string filename) const;
  void read(std::string filename);

  // to minimize the number of dependencies in ANGLE (even FFTW!),
  // commenting this out for now
  //void FFT(FIELD_2D& real, FIELD_2D& im);
  //void shiftFFT();

  // set this field to the result of convolving filter and input
  void convolve(const FIELD_2D& filter, const FIELD_2D& input);

  void resizeAndWipe(int xRes, int yRes);

  // overloaded operators
  FIELD_2D& operator=(const REAL& alpha);
  FIELD_2D& operator=(const FIELD_2D& A);
  FIELD_2D& operator*=(const REAL& alpha);
  FIELD_2D& operator/=(const REAL& alpha);
  FIELD_2D& operator+=(const REAL& alpha);
  FIELD_2D& operator-=(const FIELD_2D& input);
  FIELD_2D& operator+=(const FIELD_2D& input);
  FIELD_2D& operator*=(const FIELD_2D& input);
  FIELD_2D& operator/=(const FIELD_2D& input);

  // sum of all entries
  REAL sum();
  
  // Compute the elements of the vertical derivative convolution kernel
  void verticalDerivativeKernel(double kMax = 10, double dk = 0.01, double sigma = 1.0, double L = 0);
  
  // Compute a radial Bessel function
  void radialBessel();
  
  // set to a bessel function
  void setToBessel(float k);

  FIELD_2D nearestNeighborUpsample(int factor);

  // set to a checkboard for debugging
  void setToCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a checkboard for debugging
  void setToRampedCheckerboard(int xChecks = 10, int yChecks = 10);
  
  // set to a ramp for debugging
  void setToRampX();
  
  // set to a ramp for debugging
  void setToRampY();

  // pass a field to fieldViewer2D
  static void fieldViewer(const FIELD_2D& field, std::string name);

  // get the projection of the field in the x direction
  VECTOR projectionX();

  // return a field for the Laplacian of this field
  FIELD_2D laplacian();
  
  // return a field for the 4th order accurate Laplacian of this field
  FIELD_2D laplacian4th();
  
  // return a field for the gradient of this field
  FIELD_2D gradient();

  // return the transpose (flip x and y)
  FIELD_2D transpose() const;

  // use wavelet upsampling to double resolution
  FIELD_2D doubleRes() const;

private:
  int _xRes;
  int _yRes;
  int _totalCells;
  REAL* _data;
};

FIELD_2D operator*(const FIELD_2D& A, const REAL alpha);
FIELD_2D operator/(const FIELD_2D& A, const REAL alpha);
FIELD_2D operator+(const FIELD_2D& A, const REAL alpha);
FIELD_2D operator*(const REAL alpha, const FIELD_2D& A);
FIELD_2D operator+(const REAL alpha, const FIELD_2D& A);
FIELD_2D operator-(const FIELD_2D& A, const FIELD_2D& B);
FIELD_2D operator+(const FIELD_2D& A, const FIELD_2D& B);

} // ANGLE

#endif
