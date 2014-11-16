import numpy.ma as ma

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m

def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NAN
    """
    return ma.median(ma.masked_where(arr!=arr, arr), **kwargs)

def bottleneck_MAD(arr, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    from bottleneck import nanmedian
    import numpy as np

    if not arr.dtype.isnative:
        kind = str(arr.dtype.kind)
        sz = str(arr.dtype.itemsize)
        dt = '=' + kind + sz
        data = arr.astype(dt)
    else:
        data = arr

    if data.ndim == 1:
        d = nanmedian(data)
        m = nanmedian(ma.fabs(data - d) / c)
    else:
        d = nanmedian(data, axis=axis)
        if axis > 0:
            aswp = np.swapaxes(data,0,axis)
        else:
            aswp = data
        m = nanmedian(ma.fabs(aswp - d) / c, axis=0)

    return m
