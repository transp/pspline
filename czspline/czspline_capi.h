/**
 * $Id: czspline_capi.h 1052 2011-11-09 16:08:59Z pletzer $
 *
 *---------------------------------------------------------------------------
 * This code was developed at Tech-X (www.txcorp.com). It is free for any one
 * to use but comes with no warranty whatsoever. Use at your own risk. 
 * Thanks for reporting bugs to pletzer@txcorp.com. 
 *---------------------------------------------------------------------------
 * 
 * C callable prototypes for czspline methods
 */

#include "f77name.h"
#include "czspline_handle_size.h"
#include <stdlib.h>

#ifndef CZSPLINE_CAPI_H
#define CZSPLINE_CAPI_H

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Constructor. Allocate and set default values for spline interpolation
 * @param handle array of _ARRSZ integers
 * @param n1[, n2[, n3]] axis sizes
 * @param ier error code (0=ok)       
 */
#define _CZSPL_S F77NAME(czspline_init1)
  void _CZSPL_S(int handle[], const int *n1, 
		const int bcs1[], 
		int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_init2)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, 
		const int bcs1[], const int bcs2[], 
		int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_init3)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const int *n3,
		const int bcs1[], const int bcs2[], const int bcs3[], 
		int *ier);
#undef _CZSPL_S

/** 
 * Constructor. Allocate and set default values for hybrid interpolation
 * @param handle array of _ARRSZ integers
 * @param n1[, n2[, n3]] axis sizes
 * @param bcs1[, bcs2[, bcs3]] left/right boundary condition type 
 *        (-1=periodic, 0=not-a-knot, 1=1st derivative, 2=2nd derivative)
 *        this should be a 2-value array
 * @param ier error code (0=ok)       
 */
#define _CZSPL_S F77NAME(czshybrid_init2)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, 
		const int hspline[], 
		const int bcs1[], const int bcs2[], 
		int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czshybrid_init3)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const int *n3,
		const int hspline[], 
		const int bcs1[], const int bcs2[], const int bcs3[], 
		int *ier);
#undef _CZSPL_S

/** 
 * Constructor. Allocate and set default values for linear interpolation
 * @param handle array of _ARRSZ integers
 * @param n1[, n2[, n3]] axis sizes
 * @param bcs1[, bcs2[, bcs3]] left/right boundary condition type 
 *        (-1=periodic, 0=not-a-knot, 1=1st derivative, 2=2nd derivative)
 *        this should be a 2-value array
 * @param ier error code (0=ok)       
 */
#define _CZSPL_S F77NAME(czlinear_init1)
  void _CZSPL_S(int handle[], const int *n1, 
		int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czlinear_init2)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, 
		int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czlinear_init3)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const int *n3,
		int *ier);
#undef _CZSPL_S

/** 
 * Destructor
 * @param handle array of _ARRSZ integers
 * @param ier error code (0=ok)       
 */
#define _CZSPL_S F77NAME(czspline_free1)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_free2)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_free3)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

/** 
 * Set axes (by default axes go from (0,1) or (0, 2*pi) if periodic
 * @param handle array of _ARRSZ integers
 * @param n1[, n2[, n3]] axis sizes
 * @param x1[, x2[, x3]] axes
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_axes1)
  void _CZSPL_S(int handle[], const int *n1, const _CZSPL_P x1[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_axes2)
  void _CZSPL_S(int handle[], const int *n1, const int *n2,
		const _CZSPL_P x1[], const _CZSPL_P x2[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_axes3)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const int *n3,
		const _CZSPL_P x1[], const _CZSPL_P x2[], const _CZSPL_P x3[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Set flag for Akima interpolation (by default spline)
 * @param handle array of _ARRSZ integers
 * @param flag (1=Akima, 0=spline)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_S F77NAME(czspline_set_ishermite1)
  void _CZSPL_S(int handle[], const int *flag, int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_set_ishermite2)
  void _CZSPL_S(int handle[], const int *flag, int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_set_ishermite3)
  void _CZSPL_S(int handle[], const int *flag, int *ier);
#undef _CZSPL_S

/** 
 * Set boundary condition values
 * @param handle array of _ARRSZ integers
 * @param bcval1[, bcval2[, bcval3]] 1st or 2nd derivative values 
 *         (depending of bc type). This should be a 2-element array.
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_bcvals1)
  void _CZSPL_S(int handle[], const _CZSPL_P bcval1[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_bcvals2)
  void _CZSPL_S(int handle[], const _CZSPL_P bcval1[], const _CZSPL_P bcval2[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_set_bcvals3)
  void _CZSPL_S(int handle[], const _CZSPL_P bcval1[], const _CZSPL_P bcval2[], const _CZSPL_P bcval3[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Compute spline or Hermite coefficients
 * @param handle array of _ARRSZ integers
 * @param n1[, n2[, n3]] axis sizes
 * @param f array of original data values of size n1*n2*n3 (contiguous in n1)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_setup1)
  void _CZSPL_S(int handle[], const int *n1, const _CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_setup2)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const _CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_setup3)
  void _CZSPL_S(int handle[], const int *n1, const int *n2, const int *n3,
	  const _CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Single value (point) interpolation
 * @param handle array of _ARRSZ integers
 * @param p1[, p2[, p3]] target coordinates
 * @param f interpolated value
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp1)
  void _CZSPL_S(int handle[], const _CZSPL_P *p1, 
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp2)
  void _CZSPL_S(int handle[], const _CZSPL_P *p1, const _CZSPL_P *p2,
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp3)
  void _CZSPL_S(int handle[], const _CZSPL_P *p1, const _CZSPL_P *p2, const _CZSPL_P *p3, 
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Cloud interpolation
 * @param handle array of _ARRSZ integers
 * @param k number of points
 * @param p1[, p2[, p3]] target coordinates (each of size k)
 * @param f interpolated values (of size k)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp1_cloud)
  void _CZSPL_S(int handle[], const int *k, const _CZSPL_P p1[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp2_cloud)
  void _CZSPL_S(int handle[], const int *k, const _CZSPL_P p1[], const _CZSPL_P p2[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp3_cloud)
  void _CZSPL_S(int handle[], const int *k, const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Array interpolation
 * @param handle array of _ARRSZ integers
 * @param k1[, k2[, k3]] number of points along each axis
 * @param p1[, p2[, p3]] target coordinates (of size resp. k1, k2, k3)
 * @param f interpolated values (of size k*k2*k3, contiguous in k1)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp1_array)
  void _CZSPL_S(int handle[], const int *k1,  
		const _CZSPL_P p1[],
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp2_array)
  void _CZSPL_S(int handle[], const int *k1, const int *k2,
		const _CZSPL_P p1[], const _CZSPL_P p2[],
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_interp3_array)
  void _CZSPL_S(int handle[], const int *k1, const int *k2, const int *k3, 
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Single point derivative computation
 * @param handle array of _ARRSZ integers
 * @param i1[, i2[, i3]] order of derivative along each axis 
 *                       max order depends on interpolation method
 *                       (0 <= i1, i2, i3 <= 2 typically)
 * @param p1[, p2[, p3]] target coordinate
 * @param f return value
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative1)
  void _CZSPL_S(int handle[],
		const int *i1,
		const _CZSPL_P *p1,  
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative2)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2,
		const _CZSPL_P *p1, const _CZSPL_P *p2, 
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative3)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2, const int *i3,
		const _CZSPL_P *p1, const _CZSPL_P *p2, const _CZSPL_P *p3, 
		_CZSPL_P *f, int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Derivative computation on a cloud of points
 * @param handle array of _ARRSZ integers
 * @param i1[, i2[, i3]] order of derivative along each axis 
 *                       max order depends on interpolation method
 *                       (0 <= i1, i2, i3 <= 2 typically)
 * @param k number of points
 * @param p1[, p2[, p3]] target coordinates
 * @param f return values
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative1_cloud)
  void _CZSPL_S(int handle[],
		const int *i1,
		const int *k,
		const _CZSPL_P p1[],  
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative2_cloud)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2,
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative3_cloud)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2, const int *i3,
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Derivative computation on target rectilinear mesh
 * @param handle array of _ARRSZ integers
 * @param i1[, i2[, i3]] order of derivative along each axis 
 *                       max order depends on interpolation method
 *                       (0 <= i1, i2, i3 <= 2 typically)
 * @param k1[, k2[, k3] number of points along each axis
 * @param p1[, p2[, p3]] target coordinates
 * @param f array of size k1*k2*k3 of return values (contiguous in k1)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative1_array)
  void _CZSPL_S(int handle[],
		const int *i1,
		const int *k1,
		const _CZSPL_P p1[],  
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative2_array)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2,
		const int *k1, const int *k2,
		const _CZSPL_P p1[], const _CZSPL_P p2[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_derivative3_array)
  void _CZSPL_S(int handle[],
		const int *i1, const int *i2, const int *i3,
		const int *k1, const int *k2, const int *k3,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P f[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Single point gradient computation
 * @param handle array of _ARRSZ integers
 * @param p1[, p2[, p3]] target coordinate
 * @param df ndim return values (df/dx, df/dy, ..)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient1)
  void _CZSPL_S(int handle[],
		const _CZSPL_P *p1,
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient2)
  void _CZSPL_S(int handle[],
		const _CZSPL_P *p1, const _CZSPL_P *p2, 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient3)
  void _CZSPL_S(int handle[],
		const _CZSPL_P *p1, const _CZSPL_P *p2, const _CZSPL_P *p3, 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Gradient computation on a cloud of points
 * @param handle array of _ARRSZ integers
 * @param k number of points
 * @param p1[, p2[, p3]] target coordinates
 * @param f k*ndim return values (contiguous in k)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient1_cloud)
  void _CZSPL_S(int handle[],
		const int *k,
		const _CZSPL_P p1[],  
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient2_cloud)
  void _CZSPL_S(int handle[],
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[], 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient3_cloud)
  void _CZSPL_S(int handle[],
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Gradient computation on target rectilinear mesh
 * @param handle array of _ARRSZ integers
 * @param k1[, k2[, k3] number of points along each axis
 * @param p1[, p2[, p3]] target coordinates
 * @param df array of size k1*k2*k3*ndim of return values (contiguous in k1)
 * @param ier error code (0=ok)       
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient1_array)
  void _CZSPL_S(int handle[],
		const int *k1,
		const _CZSPL_P p1[],  
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient2_array)
  void _CZSPL_S(int handle[],
		const int *k1, const int *k2,
		const _CZSPL_P p1[], const _CZSPL_P p2[], 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_gradient3_array)
  void _CZSPL_S(int handle[],
		const int *k1, const int *k2, const int *k3,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		_CZSPL_P df[], int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Save object in file
 * @param handle array of _ARRSZ integers
 * @param filename 
 * @param ier error code (0=ok)
 */
#define _CZSPL_S F77NAME(czspline_save1)
  void _CZSPL_S(int handle[], const char *filename, int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_save2)
  void _CZSPL_S(int handle[], const char *filename, int *ier, size_t sz);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_save3)
  void _CZSPL_S(int handle[], const char *filename, int *ier, size_t sz);
#undef _CZSPL_S

/** 
 * Load object from file
 * @param handle array of _ARRSZ integers
 * @param filename 
 * @param ier error code (0=ok)
 */
#define _CZSPL_S F77NAME(czspline_load1)
  void _CZSPL_S(int handle[], const char *filename, int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_load2)
  void _CZSPL_S(int handle[], const char *filename, int *ier, size_t sz);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_load3)
  void _CZSPL_S(int handle[], const char *filename, int *ier, size_t sz);
#undef _CZSPL_S

/** 
 * Test if grid is regular, ie axes are monotonically increasing
 * @param handle array of _ARRSZ integers
 * @param ier error code (0=ok)
 */
#define _CZSPL_S F77NAME(czspline_isgridregular1)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_isgridregular2)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

#define _CZSPL_S F77NAME(czspline_isgridregular3)
  void _CZSPL_S(int handle[], int *ier);
#undef _CZSPL_S

/** 
 * Test if point coordinate is in domain
 * @param handle array of _ARRSZ integers
 * @param p1[, p2[, p3]] coordinate
 * @param ier error code (0=ok)
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain1)
  void _CZSPL_S(int handle[], 
		const _CZSPL_P *p1,
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain2)
  void _CZSPL_S(int handle[], 
		const _CZSPL_P *p1, const _CZSPL_P *p2,
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain3)
  void _CZSPL_S(int handle[], 
		const _CZSPL_P *p1, const _CZSPL_P *p2, const _CZSPL_P *p3, 
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Test if all clouds point coordinates are in domain
 * @param handle array of _ARRSZ integers
 * @param k number of points
 * @param p1[, p2[, p3]] coordinates
 * @param ier error code (0=ok)
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain1_cloud)
  void _CZSPL_S(int handle[], 
		const int *k,
		const _CZSPL_P p1[],
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain2_cloud)
  void _CZSPL_S(int handle[], 
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[],
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain3_cloud)
  void _CZSPL_S(int handle[], 
		const int *k,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

/** 
 * Test if all array point coordinates are in domain
 * @param handle array of _ARRSZ integers
 * @param k1[, k2[, k3]] number of points along each axis
 * @param p1[, p2[, p3]] coordinates along each axis
 * @param ier error code (0=ok)
 */
#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain1_array)
  void _CZSPL_S(int handle[], 
		const int *k1,
		const _CZSPL_P p1[],
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain2_array)
  void _CZSPL_S(int handle[], 
		const int *k1, const int *k2,
		const _CZSPL_P p1[], const _CZSPL_P p2[],
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S

#define _CZSPL_P double
#define _CZSPL_S F77NAME(czspline_isindomain3_array)
  void _CZSPL_S(int handle[], 
		const int *k1, const int *k2, const int *k3,
		const _CZSPL_P p1[], const _CZSPL_P p2[], const _CZSPL_P p3[], 
		int *ier);
#undef _CZSPL_P
#undef _CZSPL_S


#ifdef __cplusplus
}
#endif
  
#endif /* CZSPLINE_CAPI_H */
