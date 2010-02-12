import numpy

def correlate2d(im1,im2,psd=False,psdshift=True):
    """
    Cross-correlation of two images of arbitrary size.  Returns an image
    cropped to the size of the first input image

    Options:
    psd - if true, return fft(im1)*fft(im2[::-1,::-1]), which is the power
        spectral density
    psdshift - if true, return the shifted psd so that the DC component is in
        the center of the image

    WARNING: Normalization may be arbitrary if you use the PSD
    """
    shape1 = im1.shape
    shape2 = im2.shape
    newshape = numpy.array(shape1)+numpy.array(shape2)
    centerx, centery = newshape/2.
    quarter1x, quarter1y = centerx - shape1[0]/2.,centery - shape1[1]/2.
    quarter3x, quarter3y = centerx + shape1[0]/2.,centery + shape1[1]/2.
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
            result = shiftfftmult[ quarter1x:quarter3x, quarter1y:quarter3y ] 
            return result
        else:
            print "Warning: no fft shift - not cropped!"
            return fftmult
    else:
        rifft = (numpy.fft.ifft2( fftmult )).real
        result = rifft[ quarter1x:quarter3x, quarter1y:quarter3y ] 
        return result

