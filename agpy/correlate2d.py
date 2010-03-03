import numpy

try:
    #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    import scipy
    import scipy.fftpack
    fft2 = scipy.fftpack.fft2
    ifft2 = scipy.fftpack.ifft2
except ImportError:
    fft2 = numpy.fft.fft2
    ifft2 = numpy.fft.ifft2


def correlate2d(im1,im2,psd=False,psdshift=True,crop=True,pad=True):
    """
    Cross-correlation of two images of arbitrary size.  Returns an image
    cropped to the size of the first input image

    Options:
    psd - if true, return fft(im1)*fft(im2[::-1,::-1]), which is the power
        spectral density
    psdshift - if true, return the shifted psd so that the DC component is in
        the center of the image

    WARNING: Normalization may be arbitrary if you use the PSD

    Most of this code is duplicated from convolve.py, but I don't simply call
    convolve(im1,im2[::-1,::-1]) because the PSD that is created in an internal
    step is useful.
    """

    # NAN catching
    im1[im1!=im1] = 0
    im2[im2!=im2] = 0

    shape1 = im1.shape
    shape2 = im2.shape
    if pad:
        # find ideal size (power of 2) for fft.  Can add shapes because they are tuples
        fsize = 2**numpy.ceil(numpy.log2(numpy.max(shape1+shape2)))
        newshape = numpy.array([fsize,fsize])
    else:
        newshape = numpy.array([numpy.max([shape1[0],shape2[0]]),numpy.max([shape1[1],shape2[1]])]) #numpy.array(shape1)+numpy.array(shape2)
    centerx, centery = newshape/2.
    im1quarter1x, im1quarter1y = centerx - shape1[0]/2.,centery - shape1[1]/2.
    im1quarter3x, im1quarter3y = centerx + shape1[0]/2.,centery + shape1[1]/2.
    im2quarter1x, im2quarter1y = centerx - shape2[0]/2.,centery - shape2[1]/2.
    im2quarter3x, im2quarter3y = centerx + shape2[0]/2.,centery + shape2[1]/2.
    bigim1 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim2 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim1[im1quarter1x:im1quarter3x,im1quarter1y:im1quarter3y] = im1
    bigim2[im2quarter1x:im2quarter3x,im2quarter1y:im2quarter3y] = im2[::-1,::-1]
    imfft1 = fft2(bigim1)
    imfft2 = fft2(bigim2)
    fftmult = imfft1*imfft2

    quarter1x = numpy.min([im1quarter1x,im2quarter1x])
    quarter3x = numpy.max([im1quarter3x,im2quarter3x])
    quarter1y = numpy.min([im1quarter1y,im2quarter1y])
    quarter3y = numpy.max([im1quarter3y,im2quarter3y])
    if psd:
        if psdshift: # default on
            shiftfftmult = numpy.fft.fftshift(fftmult)
            result = shiftfftmult[ quarter1x:quarter3x, quarter1y:quarter3y ] 
            return result
        else:
            print "Warning: no fft shift - not cropped!"
            return fftmult
    else:
        rifft = numpy.fft.fftshift(ifft2( fftmult ))
        if crop: # default on
            return rifft[ quarter1x:quarter3x, quarter1y:quarter3y ]
        else: 
            return rifft

