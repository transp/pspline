/**
 *---------------------------------------------------------------------------
 * This code was developed at Tech-X (www.txcorp.com). It is free for any one
 * to use but comes with no warranty whatsoever. Use at your own risk. 
 * Thanks for reporting bugs to pletzer@txcorp.com. 
 *---------------------------------------------------------------------------
 * 
 * Test C++ bindings to ezspline
 */
#include "czspline_capi.h"

// std includes
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

/** 
 * Test 1d interpolation
 */
void test1d() {

  // opaque handle
  int handle[_ARRSZ];
 
  // error code
  int ier = 0;
  // boundary condition type
  int bcs1[2]; 
  // not-a-knot
  bcs1[0] = 0; bcs1[1] = 0;
  
  // create axes
  int n1 = 11;
  double dx = 1.0 / double(n1 - 1);
  std::vector<double> x1(n1);
  for (int i = 0; i < n1; ++i) 
    x1[i] = i * dx;
  // constructor
  F77NAME(czspline_init1)(handle, &n1, bcs1, &ier);
  assert(ier == 0);
  F77NAME(czspline_set_axes1)(handle, &n1, &x1[0], &ier);
  assert(ier == 0);

  // set original data
  std::vector<double> f(n1); 
  for (int i = 0; i < n1; ++i) 
    f[i] = x1[i]*x1[i]*x1[i];
  F77NAME(czspline_setup1)(handle, &n1, &f[0], &ier);
  assert(ier == 0);
  
  // target axes
  int m1 = 101;
  double dy = 1.0 / double(m1 - 1);
  std::vector<double> y1(m1);
  for (int i = 0; i < m1; ++i) 
    y1[i] = i * dy;

  // interpolated values
  std::vector<double> g(m1);
  
  // point interpolation
  int m = m1/2;
  double y = y1[m];
  F77NAME(czspline_interp1)(handle, &y, &g[m], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[m] - y*y*y);

  // cloud interpolation
  F77NAME(czspline_interp1_cloud)(handle, &m1, &y1[0], &g[0], &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i = 0; i < m1; ++i)
    error_cloud += std::abs(g[i] - y1[i]*y1[i]*y1[i]);
  error_cloud /= double(m1);

  // array interpolation
  F77NAME(czspline_interp1_array)(handle, &m1, &y1[0], &g[0], &ier);
  assert(ier == 0);
  double error_array = 0.0;
  for (int i = 0; i < m1; ++i)
    error_array += std::abs(g[i] - y1[i]*y1[i]*y1[i]);
  error_array /= double(m1);

  // save to file
  std::string fname = "test1d.nc";
  F77NAME(czspline_save1)(handle, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // load from file
  int handle2[_ARRSZ];
  F77NAME(czspline_load1)(handle2, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // clean up
  F77NAME(czspline_free1)(handle, &ier);
  assert(ier == 0);
  F77NAME(czspline_free1)(handle2, &ier);
  assert(ier == 0);

  std::cout << "test1d interpolation errors\n";
  std::cout << "point: " << error_point << '\n';
  std::cout << "cloud: " << error_cloud << '\n';
  std::cout << "array: " << error_array << '\n';
}

/** 
 * Test 2d interpolation
 */
void test2d() {

  // opaque handle
  int handle[_ARRSZ];

  // error code
  int ier = 0;
  // boundary condition type
  int bcs1[2]; 
  // not-a-knot
  bcs1[0] = 0; bcs1[1] = 0;
  int bcs2[2]; 
  // periodic
  bcs2[0] = -1; bcs2[1] = -1;
  
  // create axes
  int n1 = 11;
  int n2 = 12;
  double pi = 3.1415926535897931;
  double dx1 = 1.0 / double(n1 - 1);
  double dx2 = 2.0 * pi / double(n2 - 1);
  std::vector<double> x1(n1);
  for (int i = 0; i < n1; ++i) 
    x1[i] = i * dx1;
  std::vector<double> x2(n2);
  for (int i = 0; i < n2; ++i) 
    x2[i] = i * dx2;
  // constructor
  F77NAME(czspline_init2)(handle, &n1, &n2, bcs1, bcs2, &ier);
  assert(ier == 0);
  F77NAME(czspline_set_axes2)(handle, &n1, &n2, &x1[0], &x2[0], &ier);
  assert(ier == 0);

  // set original data
  std::vector<double> f(n1*n2);
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i1 = 0; i1 < n1; ++i1) {
      int i = i1 + n1 * i2;
      f[i] = x1[i1]*x1[i1]*x1[i1] * cos(x2[i2]);
    }
  }
  F77NAME(czspline_setup2)(handle, &n1, &n2, &f[0], &ier);
  assert(ier == 0);
  
  // target axes
  int m1 = 101;
  int m2 = 102;
  double dy1 = 1.0 / double(m1 - 1);
  double dy2 = 2.0 * pi / double(m2 - 1);
  std::vector<double> y1(m1);
  for (int i = 0; i < m1; ++i) 
    y1[i] = i * dy1;
  std::vector<double> y2(m2);
  for (int i = 0; i < m2; ++i) 
    y2[i] = i * dy2;

  // interpolated values
  std::vector<double> g(m1 * m2);
  
  // point interpolation
  int i1p = m1/2;
  int i2p = m2/2;
  double y1p = y1[i1p];
  double y2p = y2[i2p];
  int i = i1p + m1 * i2p;
  F77NAME(czspline_interp2)(handle, &y1p, &y2p, &g[i], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[i] - y1p*y1p*y1p * cos(y2p));

  // cloud interpolation (m1 points)
  i1p = 0;
  i = i1p + i2p * m1;
  std::vector<double> y2pp(m1, y2p);
  F77NAME(czspline_interp2_cloud)(handle, &m1, &y1[i1p], &y2pp[0], 
				  &g[i], &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i1 = 0; i1 < m1; ++i1)
    error_cloud += std::abs(g[i + i1] - y1[i1]*y1[i1]*y1[i1]*cos(y2[i2p]));
  error_cloud /= double(m1);

  // array interpolation
  F77NAME(czspline_interp2_array)(handle, &m1, &m2, 
				  &y1[0], &y2[0], &g[0], &ier);
  assert(ier == 0);
  double error_array = 0.0;
  for (int i2 = 0; i2 < m2; ++i2) {
    for (int i1 = 0; i1 < m1; ++i1) {
      int i = i1 + i2 * m1;
      error_array += std::abs(g[i] - y1[i1]*y1[i1]*y1[i1] * cos(y2[i2]));
    }
  }
  error_array /= double(m1 * m2);

  // save to file
  std::string fname = "test2d.nc";
  F77NAME(czspline_save2)(handle, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // load from file
  int handle2[_ARRSZ];
  F77NAME(czspline_load2)(handle2, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // clean up
  F77NAME(czspline_free2)(handle, &ier);
  assert(ier == 0);
  F77NAME(czspline_free2)(handle2, &ier);
  assert(ier == 0);

  std::cout << "test2d interpolation errors\n";
  std::cout << "point: " << error_point << '\n';
  std::cout << "cloud: " << error_cloud << '\n';
  std::cout << "array: " << error_array << '\n';
}

/** 
 * Test 3d interpolation
 */
void test3d() {

  // opaque handle
  int handle[_ARRSZ];

  // error code
  int ier = 0;
  // boundary condition type
  int bcs1[2]; 
  // not-a-knot
  bcs1[0] = 0; bcs1[1] = 0;
  int bcs2[2]; 
  // periodic
  bcs2[0] = -1; bcs2[1] = -1;
  int bcs3[2];
  // slope and 2nd derivative
  bcs3[0] = +1; bcs3[1] = +2;

  // create axes
  int n1 = 11;
  int n2 = 12;
  int n3 = 13;
  double pi = 3.1415926535897931;
  double dx1 = 1.0 / double(n1 - 1);
  double dx2 = 2.0 * pi / double(n2 - 1);
  double dx3 = 2.0 * pi / double(n3 - 1);
  std::vector<double> x1(n1);
  for (int i = 0; i < n1; ++i) 
    x1[i] = i * dx1;
  std::vector<double> x2(n2);
  for (int i = 0; i < n2; ++i) 
    x2[i] = i * dx2;
  std::vector<double> x3(n3);
  for (int i = 0; i < n3; ++i) 
    x3[i] = i * dx3;
  // constructor
  F77NAME(czspline_init3)(handle, &n1, &n2, &n3, 
			  bcs1, bcs2, bcs3, &ier);
  assert(ier == 0);
  // boundary conditions (values for periodic and not-a-knot
  // are not used but must still be passed to the setter)
  double bcval1[2];
  double bcval2[2];
  double bcval3[2];
  bcval3[0] = 1.0; // df/dx3 @ x3=0
  bcval3[1] = 0.0; // d^2/dx3^2 @ x3=2*pi
  F77NAME(czspline_set_bcvals3)(handle, bcval1, bcval2, bcval3, &ier);
  assert(ier == 0);
  // set axes (necessary unless (0,..1) for anything but periodic BCs
  // or (0..2*pi) for periodic BCs
  F77NAME(czspline_set_axes3)(handle, &n1, &n2, &n3, 
			      &x1[0], &x2[0], &x3[0], &ier);
  assert(ier == 0);

  // set original data
  std::vector<double> f(n1*n2*n3);
  for (int i3 = 0; i3 < n3; ++i3) {
    for (int i2 = 0; i2 < n2; ++i2) {
      for (int i1 = 0; i1 < n1; ++i1) {
	int i = i1 + n1*(i2 + n2*i3);
	f[i] = x1[i1]*x1[i1]*x1[i1] * cos(x2[i2]) * sin(x3[i3]);
      }
    }
  }
  F77NAME(czspline_setup3)(handle, &n1, &n2, &n3, 
			   &f[0], &ier);
  assert(ier == 0);
  
  // target axes
  int m1 = 101;
  int m2 = 102;
  int m3 = 103;
  double dy1 = 1.0 / double(m1 - 1);
  double dy2 = 2.0 * pi / double(m2 - 1);
  double dy3 = 2.0 * pi / double(m3 - 1);
  std::vector<double> y1(m1);
  for (int i = 0; i < m1; ++i) 
    y1[i] = i * dy1;
  std::vector<double> y2(m2);
  for (int i = 0; i < m2; ++i) 
    y2[i] = i * dy2;
  std::vector<double> y3(m3);
  for (int i = 0; i < m3; ++i) 
    y3[i] = i * dy3;

  // interpolated values
  std::vector<double> g(m1 * m2 * m3);
  
  // point interpolation
  int i1p = m1/2;
  int i2p = m2/2;
  int i3p = m3/2;
  double y1p = y1[i1p];
  double y2p = y2[i2p];
  double y3p = y3[i3p];
  int i = i1p + m1*(i2p + m2*i3p);
  F77NAME(czspline_interp3)(handle, &y1p, &y2p, &y3p, 
			    &g[i], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[i] - y1p*y1p*y1p * cos(y2p) * sin(y3p));

  // cloud interpolation (m1 points along x1 axis)
  i1p = 0;
  i = i1p + m1*(i2p + m2*i3p);
  std::vector<double> y2pp(m1, y2p);
  std::vector<double> y3pp(m1, y3p);
  F77NAME(czspline_interp3_cloud)(handle, &m1, 
				  &y1[i1p], &y2pp[0], &y3pp[0],
				  &g[i], &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i1 = 0; i1 < m1; ++i1)
    error_cloud += std::abs(g[i + i1] - 
			    y1[i1]*y1[i1]*y1[i1]*cos(y2[i2p])*sin(y3[i3p]));
  error_cloud /= double(m1);

  // array interpolation
  F77NAME(czspline_interp3_array)(handle, &m1, &m2, &m3,
				  &y1[0], &y2[0], &y3[0],
				  &g[0], &ier);
  assert(ier == 0);
  double error_array = 0.0;
  for (int i3 = 0; i3 < m3; ++i3) {
    for (int i2 = 0; i2 < m2; ++i2) {
      for (int i1 = 0; i1 < m1; ++i1) {
	int i = i1 + m1*(i2 + m2*i3);
	error_array += std::abs(g[i] - 
				y1[i1]*y1[i1]*y1[i1] * cos(y2[i2]) * sin(y3[i3]));
      }
    }
  }
  error_array /= double(m1 * m2 * m3);

  // save to file
  std::string fname = "test3d.nc";
  F77NAME(czspline_save3)(handle, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // load from file
  int handle2[_ARRSZ];
  F77NAME(czspline_load3)(handle2, fname.c_str(), &ier, fname.size());
  assert(ier == 0);

  // clean up
  F77NAME(czspline_free3)(handle, &ier);
  assert(ier == 0);
  F77NAME(czspline_free3)(handle2, &ier);
  assert(ier == 0);

  std::cout << "test3d interpolation errors\n";
  std::cout << "point: " << error_point << '\n';
  std::cout << "cloud: " << error_cloud << '\n';
  std::cout << "array: " << error_array << '\n';
}

/**
 * Main driver
 */
int main() {
    test1d();
    test2d();
    test3d();
    return 0;
}
