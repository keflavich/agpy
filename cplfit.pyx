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

from pylab import *
import numpy
cimport numpy
DTYPE=numpy.float
ctypedef numpy.float_t DTYPE_t

def plfit(x,nosmall=False,finite=False,quiet=False,silent=False):
    """
    (actually, a cython version of the python implementation....)
    A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
    from http://www.santafe.edu/~aaronc/powerlaws/

    See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
    in empirical data" SIAM Review, to appear (2009). (arXiv:0706.1062)
    http://arxiv.org/abs/0706.1062
    """
    x = asarray(x,dtype=DTYPE)
    cdef numpy.ndarray[DTYPE_t,ndim=1] xmins = unique(x)[:-1]
    #xmins = xmins[1:-1]
    cdef numpy.ndarray[DTYPE_t,ndim=1] dat = xmins * 0 
    cdef numpy.ndarray[DTYPE_t,ndim=1] av  = xmins * 0
    cdef numpy.ndarray[DTYPE_t,ndim=1] z = sort(x)
    cdef numpy.ndarray[DTYPE_t,ndim=1] cx,cf
    cdef int xm
    cdef float a,n,xmin
    for xm in xrange(len(xmins)):
        xmin = xmins[xm]
        z    = z[z>=xmin] 
        n    = float(len(z))
        # estimate alpha using direct MLE
        a    =  n / sum( log(z/xmin) )
        if nosmall:
            # 4. For continuous data, PLFIT can return erroneously large estimates of 
            #    alpha when xmin is so large that the number of obs x >= xmin is very 
            #    small. To prevent this, we can truncate the search over xmin values 
            #    before the finite-size bias becomes significant by calling PLFIT as
            if (a-1)/sqrt(n) > 0.1:
                #dat(xm:end) = [];
                dat = dat[:xm]
                xm = len(xmins)+1
                break
        # compute KS statistic
        cx   = arange(n)/float(n)  #data
        cf   = 1-(xmin/z)**a  # fitted
        av[xm] = n/sum(log(z/xmin))
        dat[xm] = max( abs(cf-cx) )
    D     = min(dat);
    xmin  = xmins[argmin(dat)]
    z     = x[x>=xmin]
    n     = len(z)
    alpha = 1 + n / sum( log(z/xmin) )
    if finite:
        alpha = alpha*(n-1.)/n+1./n
    if n < 50 and not finite and not silent:
        print '(PLFIT) Warning: finite-size bias may be present. n=%i' % n
    ks = max(abs( arange(n)/float(n) - (1-(xmin/z)**alpha) ))
    L = n*log((alpha-1)/xmin) - alpha*sum(log(z/xmin));
    #requires another map... Larr = arange(len(unique(x))) * log((av-1)/unique(x)) - av*sum
    _likelihood = L
    _av = av
    _xmin_kstest = dat
    _xmin = xmin
    _alpha= alpha
    _alphaerr = (alpha-1)/sqrt(n)
    _ks = ks
    _ngtx = n

    if not quiet:
        print "xmin: %g  n(>xmin): %i  alpha: %g +/- %g  Likelihood: %g  ks: %g" % (xmin,n,alpha,_alphaerr,L,ks)

    return xmin,alpha


def test_pl(xmin,alpha,x,ks,niter=1e3,**kwargs):
    """
    Monte-Carlo test to determine whether distribution is consistent with a power law

    Runs through niter iterations of a sample size identical to the input sample size.

    Will randomly select values from the data < xmin.  The number of values selected will
    be chosen from a uniform random distribution with p(<xmin) = n(<xmin)/n.

    Once the sample is created, is fit using above methods, then the best fit is used to
    compute a Kolmogorov-Smirnov statistic.  The KS stat distribution is compared to the 
    KS value for the fit to the actual data, and p = fraction of random ks values greater
    than the data ks value is computed.  If p<.1, the data may be inconsistent with a 
    powerlaw.  A data set of n(>xmin)>100 is required to distinguish a PL from an exponential,
    and n(>xmin)>~300 is required to distinguish a log-normal distribution from a PL.
    For more details, see figure 4.1 and section

    **WARNING** This can take a very long time to run!  Execution time scales as 
    niter * setsize

    """
    niter = int(niter)

    ntail = sum(x >= xmin)
    ntot = len(x)
    nnot = ntot-ntail              # n(<xmin)
    pnot = nnot/float(ntot)        # p(<xmin)
    nonpldata = x[x<xmin]
    nrandnot = sum( rand(ntot) < pnot ) # randomly choose how many to sample from <xmin
    nrandtail = ntot - nrandnot         # and the rest will be sampled from the powerlaw

    ksv = []
    for i in xrange(niter):
        # first, randomly sample from power law
        # with caveat!  
        nonplind = floor(rand(nrandnot)*nnot).astype('int')
        fakenonpl = nonpldata[nonplind]
        randarr = rand(nrandtail)
        fakepl = randarr**(1/(1-alpha)) * xmin 
        fakedata = concatenate([fakenonpl,fakepl])
        # second, fit to powerlaw
        TEST = plfit(fakedata,quiet=True,silent=True,nosmall=True,**kwargs)
        ksv.append(TEST._ks)
    
    ksv = array(ksv)
    p = (ksv>ks).sum() / float(niter)
    pval = p
    ks_rand = ksv

    print "p(%i) = %0.3f" % (niter,p)

    return p,ksv
