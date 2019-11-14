/**
 * Test C++ bindings to ezspline
 */

#include "czspline_capi.h"

// std includes
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>

/** 
 * Test 1d interpolation
 */
void test1d() {

  // error code
  int ier = 0;
  // boundary condition type
  int bcs1[2];
  // not-a-knot
  bcs1[0] = 0;
  bcs1[1] = 0;
  
  // create axes
  int n = 11;
  double dx = 1.0 / double(n - 1);
  double x[n];
  for (int i = 0; i < n; ++i)
    x[i] = i * dx;

  // constructor
  czspline_init1(&n, bcs1, &ier);
  assert(ier == 0);
  czspline_set_axes1(&n, x, &ier);
  assert(ier == 0);

  // set original data
  double f[n];
  for (int i = 0; i < n; ++i)
    f[i] = x[i]*x[i]*x[i];
  czspline_setup1(&n, f, &ier);
  assert(ier == 0);

  // target axes
  int m = 101;
  double dy = 1.0 / double(m - 1);
  double y[m];
  for (int i = 0; i < m; ++i)
    y[i] = i * dy;

  // interpolated values
  double g[m];

  // point interpolation
  int l = 1/2;
  double yp = y[l];
  czspline_interp1(&yp, &g[l], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[l] - yp*yp*yp);

  // cloud interpolation
  czspline_interp1_cloud(&m, y, g, &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i = 0; i < m; ++i)
    error_cloud += std::abs(g[i] - y[i]*y[i]*y[i]);
  error_cloud /= double(m);

  // array interpolation
  czspline_interp1_array(&m, y, g, &ier);
  assert(ier == 0);
  double error_array = 0.0;
  for (int i = 0; i < m; ++i)
    error_array += std::abs(g[i] - y[i]*y[i]*y[i]);
  error_array /= double(m);

#ifdef _NETCDF
  // save to file
  char fname[10] = "test1d.nc";
  czspline_save1(fname, &ier);
  assert(ier == 0);
#endif

  // clean up
  czspline_free1(&ier);
  assert(ier == 0);

#ifdef _NETCDF
  // load from file
  czspline_load1(fname, &ier);
  assert(ier == 0);

  czspline_free1(&ier);
  assert(ier == 0);
#endif

  std::cout << std::endl;
  std::cout << "1D test interpolation errors" << std::endl;
  std::cout << "point: " << error_point << std::endl;
  std::cout << "cloud: " << error_cloud << std::endl;
  std::cout << "array: " << error_array << std::endl;
}


/** 
 * Test 2d interpolation
 */
void test2d() {

  // error code
  int ier = 0;
  // boundary condition type
  int bcs1[2]; 
  // not-a-knot
  bcs1[0] = 0;
  bcs1[1] = 0;
  int bcs2[2]; 
  // periodic
  bcs2[0] = -1;
  bcs2[1] = -1;
  
  // create axes
  int n1 = 11;
  int n2 = 12;
  double pi = 3.1415926535897931;
  double dx1 = 1.0 / double(n1 - 1);
  double dx2 = 2.0 * pi / double(n2 - 1);
  double x1[n1];
  for (int i = 0; i < n1; ++i)
    x1[i] = i * dx1;
  double x2[n2];
  for (int i = 0; i < n2; ++i)
    x2[i] = i * dx2;

  // constructor
  czspline_init2(&n1, &n2, bcs1, bcs2, &ier);
  assert(ier == 0);
  czspline_set_axes2(&n1, &n2, x1, x2, &ier);
  assert(ier == 0);

  // set original data
  double f[n1*n2];
  for (int i2 = 0; i2 < n2; ++i2) {
    for (int i1 = 0; i1 < n1; ++i1) {
      int i = i1 + n1 * i2;
      f[i] = x1[i1]*x1[i1]*x1[i1] * cos(x2[i2]);
    }
  }
  czspline_setup2(&n1, &n2, f, &ier);
  assert(ier == 0);

  // target axes
  int m1 = 101;
  int m2 = 102;
  double dy1 = 1.0 / double(m1 - 1);
  double dy2 = 2.0 * pi / double(m2 - 1);
  double y1[m1];
  for (int i = 0; i < m1; ++i)
    y1[i] = i * dy1;
  double y2[m2];
  for (int i = 0; i < m2; ++i)
    y2[i] = i * dy2;

  // interpolated values
  double gm1[m1];
  double g[m1*m2];

  // point interpolation
  int i1p = m1/2;
  int i2p = m2/2;
  double y1p = y1[i1p];
  double y2p = y2[i2p];
  int i = i1p + m1 * i2p;
  czspline_interp2(&y1p, &y2p, &g[i], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[i] - y1p*y1p*y1p * cos(y2p));

  // cloud interpolation (m1 points)
  double y2pp[m1];
  for (int i = 0; i < m1; ++i)
    y2pp[i] = y2p;
  czspline_interp2_cloud(&m1, y1, y2pp, gm1, &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i1 = 0; i1 < m1; ++i1)
    error_cloud += std::abs(gm1[i1] - y1[i1]*y1[i1]*y1[i1]*cos(y2[i2p]));
  error_cloud /= double(m1);

  // array interpolation
  czspline_interp2_array(&m1, &m2, y1, y2, g, &ier);
  assert(ier == 0);
  double error_array = 0.0;
  for (int i2 = 0; i2 < m2; ++i2) {
    for (int i1 = 0; i1 < m1; ++i1) {
      int i = i1 + i2 * m1;
      error_array += std::abs(g[i] - y1[i1]*y1[i1]*y1[i1] * cos(y2[i2]));
    }
  }
  error_array /= double(m1 * m2);

#ifdef _NETCDF
  // save to file
  char fname[10] = "test2d.nc";
  czspline_save2(fname, &ier);
  assert(ier == 0);
#endif

  // clean up
  czspline_free2(&ier);
  assert(ier == 0);

#ifdef _NETCDF
  // load from file
  czspline_load2(fname, &ier);
  assert(ier == 0);

  czspline_free2(&ier);
  assert(ier == 0);
#endif

  std::cout << std::endl;
  std::cout << "2D test interpolation errors" << std::endl;
  std::cout << "point: " << error_point << std::endl;
  std::cout << "cloud: " << error_cloud << std::endl;
  std::cout << "array: " << error_array << std::endl;
}

/** 
 * Test 3d interpolation
 */
void test3d() {

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
  double x1[n1];
  for (int i = 0; i < n1; ++i) 
    x1[i] = i * dx1;
  double x2[n2];
  for (int i = 0; i < n2; ++i) 
    x2[i] = i * dx2;
  double x3[n3];
  for (int i = 0; i < n3; ++i) 
    x3[i] = i * dx3;

  // constructor
  czspline_init3(&n1, &n2, &n3, bcs1, bcs2, bcs3, &ier);
  assert(ier == 0);

  // boundary conditions (values for periodic and not-a-knot
  // are not used but must still be passed to the setter)
  double bcval1[2];
  double bcval2[2];
  double bcval3[2];
  bcval3[0] = 1.0; // df/dx3 @ x3=0
  bcval3[1] = 0.0; // d^2/dx3^2 @ x3=2*pi
  czspline_set_bcvals3(bcval1, bcval2, bcval3, &ier);
  assert(ier == 0);

  // set axes (necessary unless (0,..1) for anything but periodic BCs
  // or (0..2*pi) for periodic BCs
  czspline_set_axes3(&n1, &n2, &n3, x1, x2, x3, &ier);
  assert(ier == 0);

  // set original data
  double f[n1*n2*n3];
  for (int i3 = 0; i3 < n3; ++i3) {
    for (int i2 = 0; i2 < n2; ++i2) {
      for (int i1 = 0; i1 < n1; ++i1) {
	int i = i1 + n1*(i2 + n2*i3);
	f[i] = x1[i1]*x1[i1]*x1[i1] * cos(x2[i2]) * sin(x3[i3]);
      }
    }
  }
  czspline_setup3(&n1, &n2, &n3, f, &ier);
  assert(ier == 0);
  
  // target axes
  int m1 = 101;
  int m2 = 102;
  int m3 = 103;
  double dy1 = 1.0 / double(m1 - 1);
  double dy2 = 2.0 * pi / double(m2 - 1);
  double dy3 = 2.0 * pi / double(m3 - 1);
  double y1[m1];
  for (int i = 0; i < m1; ++i) 
    y1[i] = i * dy1;
  double y2[m2];
  for (int i = 0; i < m2; ++i) 
    y2[i] = i * dy2;
  double y3[m3];
  for (int i = 0; i < m3; ++i) 
    y3[i] = i * dy3;

  // interpolated values
  double gm1[m1];
  double g[m1 * m2 * m3];
  
  // point interpolation
  int i1p = m1/2;
  int i2p = m2/2;
  int i3p = m3/2;
  double y1p = y1[i1p];
  double y2p = y2[i2p];
  double y3p = y3[i3p];
  int i = i1p + m1*(i2p + m2*i3p);
  czspline_interp3(&y1p, &y2p, &y3p, &g[i], &ier);
  assert(ier == 0);
  double error_point = std::abs(g[i] - y1p*y1p*y1p * cos(y2p) * sin(y3p));

  // cloud interpolation (m1 points along x1 axis)
  i1p = 0;
  i = i1p + m1*(i2p + m2*i3p);
  double y2pp[m1];
  double y3pp[m1];
  for (int i = 0; i < m1; ++i) {
    y2pp[i] = y2p;
    y3pp[i] = y3p;
  }
  czspline_interp3_cloud(&m1, y1, y2pp, y3pp, gm1, &ier);
  assert(ier == 0);
  double error_cloud = 0.0;
  for (int i1 = 0; i1 < m1; ++i1)
    error_cloud += std::abs(gm1[i1] - 
 			    y1[i1]*y1[i1]*y1[i1]*cos(y2[i2p])*sin(y3[i3p]));
  error_cloud /= double(m1);

  // array interpolation
  czspline_interp3_array(&m1, &m2, &m3, y1, y2, y3, g, &ier);
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

#ifdef _NETCDF
  // save to file
  char fname[10] = "test3d.nc";
  czspline_save3(fname, &ier);
  assert(ier == 0);
#endif

  // clean up
  czspline_free3(&ier);
  assert(ier == 0);

#ifdef _NETCDF
  // load from file
  czspline_load3(fname, &ier);
  assert(ier == 0);

  czspline_free3(&ier);
  assert(ier == 0);
#endif

  std::cout << std::endl;
  std::cout << "3D test interpolation errors" << std::endl;
  std::cout << "point: " << error_point << std::endl;
  std::cout << "cloud: " << error_cloud << std::endl;
  std::cout << "array: " << error_array << std::endl;
}

/**
 * Main driver
 */
int main() {
  std::cout << std::endl << "czspline_test" << std::endl;
  test1d();
  test2d();
  test3d();
  std::cout << std::endl;
  return 0;
}
