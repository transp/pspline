CZSPLINE: C-callable interface to EZSPLINE routines
====================================================

Although the kernel of PSPLINE uses a FORTRAN 77 compatible interface, 
the preferred way to invoke PSPLINE routines is through the Fortran90 
interface routines EZSPLINE.

Unfortunately, because of lack of compatibility between Fortran90 and C,
it is not possible to use the EZSPLINE layer from another language than
Fortran. For users wanting to invoke PSPLINE/EZSPLINE from C, C++, 
Python, and other languages, we have developed a C-compatible layer 
(CZSPLINE), which preserves the object oriented flavor of the original 
EZSPLINE interface. For the most part, the EZSPLINE documentation can be 
translated to CZSPLINE is straightforward manner thanks to a nearly 
one-to-one mapping between the EZSPLINE routines and their corresponding
CZSPLINE functionality. CZSPLINE is a nearly complete layer written in 
Fortran that uses EZSPLINE underneath.

CZSPLINE has been completely rewritten to use the Fortran instrinsic
iso_c_bindings in 10/2019.
