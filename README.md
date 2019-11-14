# pspline

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/transp/pspline?include_prereleases)

## About

PSPLINE is a collection of Spline and Hermite interpolation tools for 1D, 2D, and 3D datasets on rectilinear grids developed as part of the TRANSP code.

The spline routines give full control over boundary conditions -- the user may specify "periodic", "not a knot", 1st derivative match, 2nd derivative match, or divided difference based boundary conditions on either end of each grid dimension. Hermite routines take as input the function value and derivatives at each grid point, giving back a representation of the function between grid points. Routines are provided for creating Hermite datasets, with appropriate boundary conditions applied. The 1D spline and Hermite routines are based on standard methods; the 2D and 3D spline or Hermite interpolation functions are constructed from 1D spline or Hermite interpolation functions in a straightforward manner. The splines are continuously twice differentiable in all directions across all grid cell boundaries and over the entire grid domain; Hermite functions are continuously once differentiable in all directions over the entire grid domain. For a representation of dimensionality N, an N-dimensional spline or Hermite function requires 2\*\*N\*(nx1\*nx2\*...\*nxN) memory words. There is also an "explicit spline" representation requiring 4\*\*N\*(nx1\*nx2\*...\*nxN) memory words: much more memory but somewhat faster computationally. Spline and Hermite interpolation functions are often much faster to evaluate than other representations using e.g. Fourier series or otherwise involving transcendental functions.

This version of PSPLINE includes EZspline, a Fortran-90 interface to the spline and Hermite routines. CZSPLINE, a c-callable interface to the EZspline routines, was rewritten in 10/2019 to use the Fortran instrinic iso_c_binding.


## Contact

transp_support@pppl.gov


## Contents

Makefile:      GNU make should be used
Makefile.def   Definitions for the Makefile

Source directories for:
   pspline     The original pspline F77 routines, reformatted into Fortran 90
   ezspline    A Fortran-90 interface to the spline and Hermite routines
   czspline    A C/C++ callable interface to ezpsline, using Fortran intrinsic iso_c_bindings

Test programs are included in the source directories.


## Configuration

See the sample configuration used at PPPL: pppl-bashrc


## Build

make             -- will build libpspline.a and the test programs
make debug       -- will build a debug version of libpspline.a
make shared      -- will build a shared library, libpspline.so
make clean       -- will remove the compiled objects and modules
make clobber     -- will run clean and remove the libs and test programs

The Makefile.def is set-up for compiling with the GNU, Intel, and Portland group compilers.


## Install

Installation will be put in a location specified by the PREFIX environment variable. The directory './build/' will be used if not specified. Installation done with:

make install     -- will install the code in $PREFIX
make uninstall   -- to delete the installation

You can specify the PREFIX with the make command:

make PREFIX=/path/to/location install
make PREFIX=/path/to/location uninstall


## Support

This project was supported by TRANSP development at the Princeton Plasma Physics Laboratory via `U.S. Department of Energy (DE-AC02-09CH11466)`.


