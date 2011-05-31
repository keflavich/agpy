import mpfit
import numpy as np

def powerfit(xax,data,err=None,alphaguess=-2.0,scaleguess=1.0,quiet=True):
    """
    Fit a power law (a line in log-space) to data as a function of x
    differs from 'plfit' because plfit fits a power law distribution, 
    this code simply fits a power law
    """
    
    logdata = np.log10(data)
    if err is None: err = np.ones(data.shape,dtype='float')

    def mpfitfun(data,err):
        def f(p,fjac=None): return [0,np.ravel(((np.log10(p[0])+np.log10(xax)*p[1])-data)/err)]
        return f
        
    mp = mpfit.mpfit(mpfitfun(logdata,err),xall=[scaleguess,alphaguess],quiet=quiet)
    fitp = mp.params

    return fitp,mp

def brokenpowerfit(xax, data, err=None, alphaguess1=0.0, alphaguess2=-2.0, scaleguess=1.0,
        breakpoint=None, quiet=True):
    """
    Fit a broken power law (a line in log-space) to data as a function of x
    differs from 'plfit' because plfit fits a power law distribution, 
    this code simply fits a power law

    This is a lot more intricate than the simple power law fit, since it involves
    fitting two power laws with different slopes

    Parameters:
    p[0] - scale
    p[1] - breakpoint
    p[2] - power 1 (xax < breakpoint)
    p[3] - power 2 (xax >= breakpoint)

    There are 5 parameters (NOT 4) returned because there are two scales that are *NOT*
    independent

    returns: scale1,scale2,breakpoint,alpha1,alpha2

    """
    
    logdata = np.log10(data)
    if err is None: err = np.ones(data.shape,dtype='float')

    def brokenpowerlaw(p):
        lowerhalf = (np.log10(p[0]) + np.log10(xax)*p[2]) * (xax < p[1])
        # find the location at which both functions must be equal
        scale2loc = np.argmin(np.abs(xax - p[1]))
        scale2 = np.log10(xax[scale2loc])*(p[2] - p[3]) + np.log10(p[0])
        upperhalf = (scale2 + np.log10(xax)*p[3]) * (xax >= p[1])
        # DEBUG print "scale1: %15g   scale2: %15g  xaxind: %5i  xaxval: %15g  lower: %15g upper: %15g" % (p[0],scale2,scale2loc,np.log10(xax[scale2loc]),lowerhalf[scale2loc-1],upperhalf[scale2loc])
        return lowerhalf+upperhalf


    def mpfitfun(data,err):
        def f(p,fjac=None): return [0,np.ravel((brokenpowerlaw(p)-data)/err)]
        return f

    if breakpoint is None: breakpoint = np.median(xax)

    parinfo = [{}, {'mpminstep':xax.min(),'mpmaxstep':xax.max(),'step':xax.min()}, {}, {}]
        
    mp = mpfit.mpfit(mpfitfun(logdata,err),xall=[scaleguess,breakpoint,alphaguess1,alphaguess2],quiet=quiet,parinfo=parinfo)
    fitp = mp.params

    scale2loc = np.argmin(np.abs(xax - fitp[1]))
    scale2 = 10**( np.log10(xax[scale2loc])*(fitp[2] - fitp[3]) + np.log10(fitp[0]) )
    fitp = np.array( [fitp[0],scale2] + fitp[1:].tolist() )

    return fitp,mp
