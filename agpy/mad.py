import numpy as N
#from scipy.stats import norm, median
#from scipy.stats.stats import nanmedian,_nanmedian
 	
def MAD(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    good = (a==a)
    a = N.asarray(a, N.float64)
    if a.ndim == 1:
        d = N.median(a[good])
        m = N.median(N.fabs(a[good] - d) / c)
    else:
        d = N.median(a[good], axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(a[good],0,axis)
        else:
            aswp = a[good]
        m = N.median(N.fabs(aswp - d) / c, axis=0)

    return m

def nanmedian(arr):
    """
    Returns median ignoring NAN
    """
    return N.median(arr[arr==arr])
