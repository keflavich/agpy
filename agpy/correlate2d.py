import numpy

def correlate2d(im1,im2,psd=False,psdshift=True):
    """
    Cross-correlation of two images of arbitrary size

    Options:
    psd - if true, return fft(im1)*fft(im2[::-1,::-1]), which is the power spectral density
    psdshift - if true, return the shifted psd so that the DC component is in the center of the image

    WARNING: Normalization may be arbitrary if you use the PSD
    """
    shape1 = im1.shape
    shape2 = im2.shape
    newshape = numpy.array(shape1)+numpy.array(shape2)
    centerx, centery = newshape/2.
    quarterx, quartery = newshape/4.
    bigim1 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim2 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim1[:im1.shape[0],:im1.shape[1]] = im1
    # cross-correlation instead of convolution: reverse the array
    bigim2[:im2.shape[0],:im2.shape[1]] = im2[::-1,::-1] 
    fft1 = numpy.fft.fft2(bigim1)
    fft2 = numpy.fft.fft2(bigim2)
    fftmult = fft1*fft2
    if psd:
        if psdshift:
            shiftfftmult = numpy.fft.fftshift(fftmult)
            result = shiftfftmult[ quarterx:centerx+quarterx, quartery:centery+quartery ] 
            return result
        else:
            print "Warning: no fft shift - not cropped!"
            return fftmult
    else:
        rifft = (numpy.fft.ifft2( fftmult )).real
        result = rifft[ quarterx:centerx+quarterx, quartery:centery+quartery ] 
        return result

