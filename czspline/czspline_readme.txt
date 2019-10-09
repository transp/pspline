CZSPLINE: C-callable interface to EZSPLINE routines
====================================================

Although the kernel of PSPLINE uses a FORTRAN 77 compatible interface, the 
preferred way to invoke PSPLINE routines is through the Fortran90 interface
routines EZSPLINE, which have an object oriented flavor and are heavily 
overloaded both in their functionality and by supporting real*4/real*8 
precision.

Unfortunately, because of lack of compatibility between Fortran90 and C, it
is not possible to use the EZSPLINE layer from another language than Fortran. 
For users wanting to invoke PSPLINE/EZSPLINE from C, C++, Python, and other languages, 
we have developed a C-compatible layer (CZSPLINE), which preserves the object
oriented flavor of the original EZSPLINE interface. For the most part, 
the EZSPLINE documentation can be translated to CZSPLINE is straightforward manner 
thanks to a nearly one-to-one mapping between the EZSPLINE routines and their 
corresponding CZSPLINE functionality. 

CZSPLINE is a nearly complete layer written in Fortran that uses 
EZSPLINE underneath. Keep in mind, however, that since FORTRAN 77/C don't 
support method overloading, the name of the CZSPLINE routines are those of
the EZSPLINE routines but extended to allow the C compiler to distinguish between
different precisions (r4 or r8) and functionality (eg, whether the interpolation
acts on 1, 2, or 3 dimensional data). 

EZSPLINE routines expect a derived type argument (the spline object). In 
CZSPLINE, this object maps to an array of integers using the pointer to 
opaque container technique [Pletzer08]. 

Comments, requests, and bug reports can be addressed to pletzer@txcorp.com.
 
Rules to map an ezspline routine to a czspline routine
------------------------------------------------------

All CZSPLINE methods use lower case names in C.

C <--> FORTRAN 77 name mangling is taken care (no need to add an underscore).

All methods are subroutines. Functions in EZSPLINE are void procedures in 
CZSPLINE, their return value being the last argument before the error code.

The following rule applies to translate an EZSPLINE routine to a CZSPLINE
corresponding routine

EZSPLINE            CZSPLINE

ezspline_<METHOD>   czspline_<METHOD><DIM>_[<MTYPE>_]<PRECISION>

where <METHOD> is one of init, free, setup, interp, derivative, gradient....,
<DIM> is either 1, 2, or 3 and refers to the data dimensionality, and 
<PRECISION> is either r4 or r8. Some routines have an additional descriptor
<MTYPE> when <DIM> and <PRECISION> are insufficient to resolve the routine
(eg czspline_interp2_cloud_r8 vs czspline_interp2_arry_r8)

Argument list
-------------

In C all arguments must be provided in pre-established order. The following 
convention applies:

First argument: int handle[] opaque handle (list of integers). At least 12 
integers are required.

Example: 
EZSPLINE: type(ezspline2) :: spl
CZSPLINE: int spl[12];

Last argument: int *ier error code (0=OK)

All other arguments: are routine dependent primitive types passed by reference

Const'ness is preseved, ie intent(in) arguments on the Fortran side have the
 "const" on the C side.

When passing arrays, the size must be provided in CZSPLINE except when these
are constant. The array sizes precede the array in the argument list.

Example:
EZSPLINE: call ezspline_interp(spl, n, p1, p1, f, ier)
CZSPLINE: czspline_interp2_cloud_r8(int spl[], const int *n, 
	    const double p1[], const double p2[], double f[],
            int *ier);

Additional setter methods
-------------------------

Data members of spline objects cannot be set directly in C. We are providing 
convenience setter methods for that purpose. These have the name 

czspline_set_<MEMBER>_<PRECISION>

Example: czspline_set_ishermite1_r4(int spl[], const int *flag, int *ier);

Axes should be set using 

czspline_set_axes<DIM>_<PRECISION>

which has the same functionality as spl%x1=..., spl%x2=...

List of CZSPLINE procedures
---------------------------

grep F77NAME czspline_capi.h | perl -ne 's/#define _S F77NAME\(//;s/\)//;print;'

czspline_init1_r4
czspline_init2_r4
czspline_init3_r4
czspline_init1_r8
czspline_init2_r8
czspline_init3_r8
czspline_free1_r4
czspline_free2_r4
czspline_free3_r4
czspline_free1_r8
czspline_free2_r8
czspline_free3_r8
czspline_set_axes1_r4
czspline_set_axes2_r4
czspline_set_axes3_r4
czspline_set_axes1_r8
czspline_set_axes2_r8
czspline_set_axes3_r8
czspline_set_ishermite1_r4
czspline_set_ishermite2_r4
czspline_set_ishermite3_r4
czspline_set_ishermite1_r8
czspline_set_ishermite2_r8
czspline_set_ishermite3_r8
czspline_set_bcvals1_r4
czspline_set_bcvals2_r4
czspline_set_bcvals3_r4
czspline_set_bcvals1_r8
czspline_set_bcvals2_r8
czspline_set_bcvals3_r8
czspline_setup1_r4
czspline_setup2_r4
czspline_setup3_r4
czspline_setup1_r8
czspline_setup2_r8
czspline_setup3_r8
czspline_interp1_r4
czspline_interp2_r4
czspline_interp3_r4
czspline_interp1_r8
czspline_interp2_r8
czspline_interp3_r8
czspline_interp1_cloud_r4
czspline_interp2_cloud_r4
czspline_interp3_cloud_r4
czspline_interp1_cloud_r8
czspline_interp2_cloud_r8
czspline_interp3_cloud_r8
czspline_interp1_array_r4
czspline_interp2_array_r4
czspline_interp3_array_r4
czspline_interp1_array_r8
czspline_interp2_array_r8
czspline_interp3_array_r8
czspline_derivative1_r4
czspline_derivative2_r4
czspline_derivative3_r4
czspline_derivative1_r8
czspline_derivative2_r8
czspline_derivative3_r8
czspline_derivative1_cloud_r4
czspline_derivative2_cloud_r4
czspline_derivative3_cloud_r4
czspline_derivative1_cloud_r8
czspline_derivative2_cloud_r8
czspline_derivative3_cloud_r8
czspline_derivative1_array_r4
czspline_derivative2_array_r4
czspline_derivative3_array_r4
czspline_derivative1_array_r8
czspline_derivative2_array_r8
czspline_derivative3_array_r8
czspline_gradient1_r4
czspline_gradient2_r4
czspline_gradient3_r4
czspline_gradient1_r8
czspline_gradient2_r8
czspline_gradient3_r8
czspline_gradient1_cloud_r4
czspline_gradient2_cloud_r4
czspline_gradient3_cloud_r4
czspline_gradient1_cloud_r8
czspline_gradient2_cloud_r8
czspline_gradient3_cloud_r8
czspline_gradient1_array_r4
czspline_gradient2_array_r4
czspline_gradient3_array_r4
czspline_gradient1_array_r8
czspline_gradient2_array_r8
czspline_gradient3_array_r8
czspline_save1_r4
czspline_save2_r4
czspline_save3_r4
czspline_save1_r8
czspline_save2_r8
czspline_save3_r8
czspline_load1_r4
czspline_load2_r4
czspline_load3_r4
czspline_load1_r8
czspline_load2_r8
czspline_load3_r8
czspline_isgridregular1_r4
czspline_isgridregular2_r4
czspline_isgridregular3_r4
czspline_isgridregular1_r8
czspline_isgridregular2_r8
czspline_isgridregular3_r8
czspline_isindomain1_r4
czspline_isindomain2_r4
czspline_isindomain3_r4
czspline_isindomain1_r8
czspline_isindomain2_r8
czspline_isindomain3_r8
czspline_isindomain1_cloud_r4
czspline_isindomain2_cloud_r4
czspline_isindomain3_cloud_r4
czspline_isindomain1_cloud_r8
czspline_isindomain2_cloud_r8
czspline_isindomain3_cloud_r8
czspline_isindomain1_array_r4
czspline_isindomain2_array_r4
czspline_isindomain3_array_r4
czspline_isindomain1_array_r8
czspline_isindomain2_array_r8
czspline_isindomain3_array_r8

References
----------

[Pletzer08] A Pletzer, D McCune, S Muszala, S Vadlamani, S Kruger: "Exposing Fortran derived types to C and other languages", Computing in Science and Engineering, July/August 2008, p.86
