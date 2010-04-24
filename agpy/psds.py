import numpy
from correlate2d import correlate2d
from radialprofile import azimuthalAverage

try:
    #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    import scipy
    import scipy.fftpack
    fft2 = scipy.fftpack.fft2
    ifft2 = scipy.fftpack.ifft2
except ImportError:
    fft2 = numpy.fft.fft2
    ifft2 = numpy.fft.ifft2


def PSD2(image,image2=None,oned=True,return_index=True,wavenumber=False,pad=False):
    """
    Two-dimensional PSD
    oned - return radial profile of 2D PSD (i.e. mean power as a function of spatial frequency)
           freq,zz = PSD2(image); plot(freq,zz) is a power spectrum
    return_index - if true, the first return item will be the indexes
    wavenumber - if one dimensional and return_index set, will return a normalized wavenumber instead
    """
    

    image[image!=image] = 0
    if image2 is None:
        image2 = image
    psd2 = numpy.abs( correlate2d(image,image2,return_fft=True,pad=pad) ) 
    # normalization is approximately (numpy.abs(image).sum()*numpy.abs(image2).sum())

    if oned:
        #freq = 1 + numpy.arange( numpy.floor( numpy.sqrt((image.shape[0]/2)**2+(image.shape[1]/2)**2) ) ) 

        freq,zz = azimuthalAverage(psd2,returnradii=True)
        freq = freq.astype('float') + 1.0

        if return_index:
            if wavenumber:
                return len(freq)/freq,zz
            else:
                return freq/len(freq),zz
        else:
            return zz
    else:
        return psd2
