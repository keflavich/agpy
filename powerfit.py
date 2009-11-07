import mpfit
import numpy as np

def powerfit(xax,data,alphaguess=-2.0,scaleguess=1.0,quiet=True):
    """
    Fit a power law (a line in log-space) to data as a function of x
    differs from 'plfit' because plfit fits a power law distribution, 
    this code simply fits a power law
    """
    
    logdata = np.log10(data)

    def mpfitfun(data):
        def f(p,fjac=None): return [0,np.ravel((np.log10(p[0])+np.log10(xax)*p[1])-data)]
        return f
        
    mp = mpfit.mpfit(mpfitfun(logdata),xall=[scaleguess,alphaguess],quiet=quiet)
    fitp = mp.params

    return fitp,mp
