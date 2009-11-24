import numpy as N
from scipy.stats import norm, median
from scipy.stats.stats import nanmedian,_nanmedian
 	
def MAD(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    """

    a = N.asarray(a, N.float64)
    if a.ndim == 1:
        d = _nanmedian(a)
        m = _nanmedian(N.fabs(a - d) / c)
    else:
        d = nanmedian(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(0,axis)
        else:
            aswp = a
        m = nanmedian(N.fabs(aswp - d) / c, axis=0)

    return m

