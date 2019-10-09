#!/usr/bin/env python

"""
---------------------------------------------------------------------------
 This code was developed at Tech-X (www.txcorp.com). It is free for any one
 to use but comes with no warranty whatsoever. Use at your own risk. 
 Thanks for reporting bugs to pletzer@txcorp.com. 
---------------------------------------------------------------------------

czspline unit test
$Id: czspline_test.py,v 1.1 2008-05-22 16:30:24 Alex_Pletzer Exp $
"""

import numpy
from ctypes import POINTER, c_int, c_double, c_float, byref

def test1d_r8():

    """
    Test 1d interpolation
    """

    # access shared object
    tr = ctypes.CDLL(libdir + "/libTranspPhage.so", \
                     mode=ctypes.RTLD_GLOBAL)
    # opaque handle
    holder = numpy.array([0]*12)
    handle = handle.ctype.data_as(ctypes.POINTER(ctypes.c_double))
    bcs1   = numpy.array([0]*2) # not-a-knot
    ier    = ctypes.int()
    
    # construct and set axes
    x1 = numpy.arange(0.0, 1.00001, 0.1)
    tr.czspline_init1_r8_(handle,
                          byref(c_int(len(x1))),
                          bcs1.ctypes.data_as(POINTER(c_int)),
                          byref(ier))
    tr.czspline_set_axes1_r8(handle,
                             byref(c_int(len(x1))),
                             x1.ctypes.data_as(POINTER(c_double)),
                             byref(ier))
    
    # original data
    f  = x1**3

    # compute spline coefficients
    tr.czspline_setup1_r8_(handle,
                      byref(c_int(f.shape[0])),
                      f.ctypes.data_as(POINTER(c_double)),
                      byref(ier))

    # target axes
    y1 = numpy.arange(0.0, 1.00001, 0.01)
    g  = 0 * y1
    
    # array interpolation
    tr.czspline_interp1_array_r8_(handle,
                       byref(c_int(y1.shape[0])),
                                  g.ctypes.data_as(POINTER(c_double)),
                                  byref(ier))
    error_array = numpy.sum( numpy.abs( g - y1**3 ) )/len(g)

    # cloud interpolation
    tr.czspline_interp1_cloud_r8_(handle,
                       byref(c_int(y1.shape[0])),
                                  g.ctypes.data_as(POINTER(c_double)),
                                  byref(ier))
    error_cloud = numpy.sum( numpy.abs( g - y1**3 ) )/len(g)
    
    # point interpolation
    g = c_double()
    i  = len(y1)/2
    tr.czspline_interp1_r8_(handle,
                            byref(c_double(y1[i])),
                            byref(g),
                            byref(ier))
    error_point = g - y1[i]**3

    # clean up
    tr.czspline_free1_r8(handle,
                         byref(ier))
                            
                                  
    print """
test1d_r8 interpolation errors
point: %g
cloud: %g
array: %g
    """ % (error_point, error_cloud, error_array)
    
    

def test2d_r8():
    pass

def test3d_r8():
    pass

##############################################################################
def main():
    import sys
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('-l', '--libdir', action='store', type='string',
                      dest='libdir', help='location of libTranspPhage',
                      default='/usr/local/lib')
    options, args = parser.parse_args()
    libdir = options.libdir
    
    test1d_r8(libdir)
    test2d_r8(libdir)
    test3d_r8(libdir)
if __name__ == '__main__': main()
