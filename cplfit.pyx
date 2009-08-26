# intended to implement a power-law fitting routine as specified in.....
# http://www.santafe.edu/~aaronc/powerlaws/
#
# The MLE for the power-law alpha is very easy to derive given knowledge
# of the lowest value at which a power law holds, but that point is 
# difficult to derive and must be acquired iteratively.

"""
plfit.py - a python power-law fitter based on code by Aaron Clauset
http://www.santafe.edu/~aaronc/powerlaws/
http://arxiv.org/abs/0706.1062 "Power-law distributions in empirical data" 
Requires pylab (matplotlib), which requires numpy

example use:
from plfit import plfit

MyPL = plfit(mydata)
MyPL.plotpdf(log=True)


"""

#from pylab import *
cimport c_numpy as cnp
cnp.import_array()
import numpy
cimport numpy
import cmath
import time
DTYPE=numpy.float
ctypedef numpy.float_t DTYPE_t

cdef extern from "math.h":
    float log(float theta)
    float sqrt(float theta)
    float abs(float theta)
    float pow(float x,float y)

def plfit_loop(z,xmins,nosmall=True):
    """
    The internal loop of plfit.  Useful because it can
    be wrapped into plfit.py for direct comparison
    with the fortran version

    z must be sorted
    """
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] dat = xmins * 0 
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] av  = xmins * 0
    cdef int lxm = len(xmins)
    cdef float cx,cf,val,a,xmin
    cdef int xm,n,i
    n = float(len(z))
    xmin = 0
    for xm from 0 <= xm < lxm:
        z=z[1:]
        n-=1
        if xmins[xm] == xmin:
            continue
        xmin = xmins[xm]
        # estimate alpha using direct MLE
        a = 0
        for i from 0 <= i < n:
            a +=  n / ( log(z[i]/xmin) )
        if nosmall:
            # 4. For continuous data, PLFIT can return erroneously large estimates of 
            #    alpha when xmin is so large that the number of obs x >= xmin is very 
            #    small. To prevent this, we can truncate the search over xmin values 
            #    before the finite-size bias becomes significant by calling PLFIT as
            if (a-1)/sqrt(n) > 0.1:
                #dat(xm:end) = [];
                dat = dat[:xm] # this should only be called once
                break
        av[xm] = a
        dat[xm] = 0
#        dat[xm] = numpy.max(numpy.abs(cf-cx))
        for i from 0 <= i < n:
            # compute KS statistic
            cx   = <float>(i)/n  #data; cheap float conversion?
            cf   = 1.0-pow((xmin/z[i]),a)  # fitted
            val = abs(cf-cx)
            if dat[xm] < val:
                dat[xm] = val

def plfit(x,nosmall=True,finite=False,quiet=False,silent=False):
    """
    (actually, a cython version of the python implementation....)
    A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
    from http://www.santafe.edu/~aaronc/powerlaws/

    See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
    in empirical data" SIAM Review, to appear (2009). (arXiv:0706.1062)
    http://arxiv.org/abs/0706.1062
    """
    x = numpy.asarray(x,dtype=DTYPE)
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] xmins = numpy.unique(x)[1:-1]
    #xmins = xmins[1:-1]
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] dat = xmins * 0 
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] av  = xmins * 0
    cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] z = numpy.sort(x)
    #cdef cnp.ndarray[DTYPE_t,ndim=1,mode='c'] cx,cf
    cdef float cx,cf,val
    cdef int xm,n,i
    cdef int lxm = len(xmins)
    cdef float a,xmin
    #for xm in xrange(len(xmins)):
    t=time.time()
    dat,av = plfit_loop(z,xmins,nosmall=nosmall)
    if not quiet: print "CYTHON plfit executed in %f seconds" % (time.time()-t)
    xmin  = xmins[numpy.argmin(dat)]
    z     = x[x>=xmin]
    n     = len(z)
    alpha = 1 + float(n) / numpy.sum( numpy.log(z/xmin) )
    if finite:
        alpha = alpha*(n-1.)/n+1./n
    if n < 50 and not finite and not silent:
        print '(PLFIT) Warning: finite-size bias may be present. n=%i' % n
    ks = numpy.max(numpy.abs( numpy.arange(n)/float(n) - (1-(xmin/z)**alpha) ))
    L = n*numpy.log((alpha-1)/xmin) - alpha*numpy.sum(numpy.log(z/xmin));
    #requires another map... Larr = arange(len(unique(x))) * log((av-1)/unique(x)) - av*sum
    _likelihood = L
    _av = av
    _xmin_kstest = dat
    _xmin = xmin
    _alpha= alpha
    _alphaerr = (alpha-1)/numpy.sqrt(n)
    _ks = ks
    _ngtx = n

    if not quiet:
        print "xmin: %g  n(>xmin): %i  alpha: %g +/- %g  Likelihood: %g  ks: %g" % (xmin,n,alpha,_alphaerr,L,ks)

    return xmin,alpha
