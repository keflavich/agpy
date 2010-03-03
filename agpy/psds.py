import numpy
from correlate2d import correlate2d

try:
    #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    import scipy
    import scipy.fftpack
    fft2 = scipy.fftpack.fft2
    ifft2 = scipy.fftpack.ifft2
except ImportError:
    fft2 = numpy.fft.fft2
    ifft2 = numpy.fft.ifft2


def PSD2(image,image2=None,oned=True,return_index=True,wavenumber=False):
    """
    Two-dimensional PSD
    oned - return radial profile of 2D PSD (i.e. mean power as a function of spatial frequency)
           plot(freq,zz) is a power spectrum
    return_index - if true, the first return item will be the indexes
    wavenumber - if one dimensional and return_index set, will return a normalized wavenumber instead
    """
    

    image[image!=image] = 0
    if image2 is None:
        image2 = image
    psd2 = numpy.abs( correlate2d(image,image2,return_fft=True) ) / (numpy.abs(image)*numpy.abs(image2)).sum()

    xx,yy = numpy.indices(image.shape)
    rr = numpy.sqrt((xx-image.shape[0] / 2)**2+(yy-image.shape[1] / 2)**2)

    if oned:
        freq = numpy.arange( numpy.floor( numpy.sqrt((image.shape[0]/2)**2+(image.shape[1]/2)**2) ) ) 

        zz = numpy.array([ psd2[xx[rr.round()==ii],yy[rr.round()==ii]].mean() for ii in freq])

        if return_index:
            if wavenumber:
                return len(freq)/freq,zz
            else:
                return freq/len(freq),zz
        else:
            return zz
    else:
        if return_index:
            return rr,psd2
        else:
            return psd2
