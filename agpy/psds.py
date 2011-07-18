import numpy
from correlate2d import correlate2d
from radialprofile import azimuthalAverageBins,radialAverageBins

try:
    #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    import scipy
    import scipy.fftpack
    fft2 = scipy.fftpack.fft2
    ifft2 = scipy.fftpack.ifft2
except ImportError:
    fft2 = numpy.fft.fft2
    ifft2 = numpy.fft.ifft2


def PSD2(image, image2=None, oned=True, return_index=True, wavenumber=False,
        fft_pad=False, return_stddev=False, real=False, imag=False,
        binsize=1.0, radbins=1, azbins=1, radial=False, **kwargs):
    """
    Two-dimensional PSD
    oned - return radial profile of 2D PSD (i.e. mean power as a function of spatial frequency)
           freq,zz = PSD2(image); plot(freq,zz) is a power spectrum
    return_index - if true, the first return item will be the indexes
    wavenumber - if one dimensional and return_index set, will return a normalized wavenumber instead
    real - Only compute the real part of the PSD
    complex - Only compute the complex part of the PSD
    """
    

    image[image!=image] = 0
    if image2 is None:
        image2 = image
    if real:
        psd2 = numpy.real( correlate2d(image,image2,return_fft=True,fft_pad=fft_pad) ) 
    elif imag:
        psd2 = numpy.imag( correlate2d(image,image2,return_fft=True,fft_pad=fft_pad) ) 
    else:
        psd2 = numpy.abs( correlate2d(image,image2,return_fft=True,fft_pad=fft_pad) ) 
    # normalization is approximately (numpy.abs(image).sum()*numpy.abs(image2).sum())

    if radial:
        azbins,az,zz = radialAverageBins(psd2,radbins=radbins, interpnan=True, binsize=binsize, **kwargs)
        if len(zz) == 1:
            return az,zz[0]
        else:
            return az,zz

    if oned:
        #freq = 1 + numpy.arange( numpy.floor( numpy.sqrt((image.shape[0]/2)**2+(image.shape[1]/2)**2) ) ) 

        azbins,(freq,zz) = azimuthalAverageBins(psd2,azbins=azbins,interpnan=True, binsize=binsize, **kwargs)
        if len(zz) == 1: zz=zz[0]
        # the "Frequency" is the spatial frequency f = 1/x for the standard numpy fft, which follows the convention
        # A_k =  \sum_{m=0}^{n-1} a_m \exp\left\{-2\pi i{mk \over n}\right\}
        # or 
        # F_f = Sum( a_m e^(-2 pi i f x_m)  over the range m,m_max where a_m are the values of the pixels, x_m are the
        # indices of the pixels, and f is the spatial frequency
        freq = freq.astype('float')  # there was a +1.0 here before, presumably to deal with div-by-0, but that shouldn't happen and shouldn't have been "accounted for" anyway

        if return_index:
            if wavenumber:
                return_vals = list((len(freq)/freq,zz))
            else:
                return_vals = list((freq/len(freq),zz))
        else:
            return_vals = list(zz)
        if return_stddev:
            zzstd = azimuthalAverageBins(psd2,azbins=azbins,stddev=True,interpnan=True, binsize=binsize, **kwargs)
            return_vals.append(zzstd)

        return return_vals

    # else...
    return psd2
