import numpy

def downsample(myarr,factor):
    """
    Downsample a 2D array by averaging over *factor* pixels in each axis.
    Crops upper edge if the shape is not a multiple of factor.

    This code is pure numpy and should be fast.
    """
    xs,ys = myarr.shape
    crarr = myarr[:xs-(xs % int(factor)),:ys-(ys % int(factor))]
    dsarr = numpy.concatenate([[crarr[i::factor,j::factor] 
        for i in range(factor)] 
        for j in range(factor)]).mean(axis=0)
    return dsarr
