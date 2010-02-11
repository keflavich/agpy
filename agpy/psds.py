import numpy
from correlate2d import correlate2d

def PSD2(image,oned=True,return_index=True,wavenumber=False):
    """
    Two-dimensional PSD
    oned - return radial profile of 2D PSD (i.e. mean power as a function of spatial frequency)
           plot(freq,zz) is a power spectrum
    return_index - if true, the first return item will be the indexes
    wavenumber - if one dimensional and return_index set, will return a normalized wavenumber instead
    """

    image[image!=image] = 0
    #acorr = scipy.stsci.convolve.correlate2d(image,image,fft=True,mode='constant')
    acorr = correlate2d(image,image)
    psd2 = numpy.abs( numpy.fft.fftshift( numpy.fft.fft2(acorr) ) )
    #psd2 = numpy.abs( ( correlate2d(image,image,psd=True) ) )

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
