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

def convolve(im1,im2):
    """
    Convolve two images.  Returns the larger image size.  Assumes image & kernel are centered

    Options:
    """
    shape1 = im1.shape
    shape2 = im2.shape
    # find ideal size (power of 2) for fft.  Can add shapes because they are tuples
    fsize = 2**numpy.ceil(numpy.log2(numpy.max(shape1+shape2)))
    newshape = numpy.array([fsize,fsize])
    centerx, centery = newshape/2.
    im1quarter1x, im1quarter1y = centerx - shape1[0]/2.,centery - shape1[1]/2.
    im1quarter3x, im1quarter3y = centerx + shape1[0]/2.,centery + shape1[1]/2.
    im2quarter1x, im2quarter1y = centerx - shape2[0]/2.,centery - shape2[1]/2.
    im2quarter3x, im2quarter3y = centerx + shape2[0]/2.,centery + shape2[1]/2.
    bigim1 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim2 = numpy.zeros(newshape,dtype=numpy.float64)
    bigim1[im1quarter1x:im1quarter3x,im1quarter1y:im1quarter3y] = im1
    bigim2[im2quarter1x:im2quarter3x,im2quarter1y:im2quarter3y] = im2 
    imfft1 = fft2(bigim1)
    imfft2 = fft2(bigim2)
    fftmult = imfft1*imfft2

    rifft = numpy.fft.fftshift( (ifft2( fftmult )).real )
    quarter1x = numpy.min([im1quarter1x,im2quarter1x])
    quarter3x = numpy.max([im1quarter3x,im2quarter3x])
    quarter1y = numpy.min([im1quarter1y,im2quarter1y])
    quarter3y = numpy.max([im1quarter3y,im2quarter3y])
    result = rifft[ quarter1x:quarter3x, quarter1y:quarter3y ] 
    return result

def smooth(image,kernelwidth=3,kerneltype='gaussian',trapslope=None):
    """
    Returns a smoothed image using a gaussian, boxcar, or tophat kernel

    Options: 
    kerneltype - gaussian, boxcar, or tophat.  
        For a gaussian, uses a gaussian with sigma = kernelwidth (in pixels) out to 9-sigma
        A boxcar is a kernelwidth x kernelwidth square 
        A tophat is a flat circle with radius = kernelwidth
    """

    if kerneltype == 'gaussian':
        gp = 9 # gaussian precision in n_sigma
        xx,yy = numpy.indices([numpy.ceil(kernelwidth*gp),numpy.ceil(kernelwidth*gp)])
        rr = numpy.sqrt((xx-kernelwidth*gp/2.)**2+(yy-kernelwidth*gp/2.)**2)
        kernel = numpy.exp(-(rr**2)/(2*kernelwidth**2)) / (kernelwidth**2 * (2*numpy.pi))
#        if kernelwidth != numpy.round(kernelwidth):
#            print "Rounding kernel width to %i pixels" % numpy.round(kernelwidth)
#            kernelwidth = numpy.round(kernelwidth)

    if kerneltype == 'boxcar':
        print "Using boxcar kernel size %i" % numpy.ceil(kernelwidth)
        kernel = numpy.ones([numpy.ceil(kernelwidth),numpy.ceil(kernelwidth)],dtype='float64') / kernelwidth**2
    elif kerneltype == 'tophat':
        kernel = numpy.zeros(image.shape,dtype='float64')
        xx,yy = numpy.indices(image.shape)
        rr = numpy.sqrt((xx-image.shape[0]/2.)**2+(yy-image.shape[1]/2.)**2)
        kernel[rr<kernelwidth] = 1.0
        # normalize
        kernel /= kernel.sum()
    elif kerneltype == 'brickwall':
        print "Smoothing with a %f pixel brick wall filter" % kernelwidth
        xx,yy = numpy.indices(image.shape)
        rr = numpy.sqrt((xx-image.shape[0]/2.)**2+(yy-image.shape[1]/2.)**2)
        kernel = numpy.sinc(rr/kernelwidth) 
        kernel /= kernel.sum()
    elif kerneltype == 'trapezoid':
        if trapslope:
            xx,yy = numpy.indices(image.shape)
            rr = numpy.sqrt((xx-image.shape[0]/2.)**2+(yy-image.shape[1]/2.)**2)
            zz = rr.max()-(rr*trapslope)
            zz[zz<0] = 0
            zz[rr<kernelwidth] = 1.0
            kernel = zz/zz.sum()
        else:
            print "trapezoid function requires a slope"

    bad = (image != image)
    temp = image.copy()
    temp[bad] = 0
    temp = convolve(temp,kernel)
    temp[bad] = image[bad]

    return temp

