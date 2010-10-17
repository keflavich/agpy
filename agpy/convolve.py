import numpy as np

try:
    #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    import scipy
    import scipy.fftpack
    fft2 = scipy.fftpack.fft2
    ifft2 = scipy.fftpack.ifft2
    from scipy.special import j1
except ImportError:
    fft2 = np.fft.fft2
    ifft2 = np.fft.ifft2

def convolve(im1, im2, crop=True, return_fft=False, fftshift=True,
        kernel_number=1, fft_pad=True, psf_pad=False, ignore_nan=False,
        ignore_zeros=False, min_wt=1e-8):
    """
    Convolve two images.  Returns the larger image size.  Assumes image &
    kernel are centered

    *NOTE* Order matters when padding; the kernel should either be first or
    kernel_number should be set to 2.  

    Options:
    fft_pad - Default on.  Zero-pad image to the nearest 2^n
    psf_pad - Default off.  Zero-pad image to be at least the sum of the image
        sizes (in order to avoid edge-wrapping when smoothing)
    crop - Default on.  Return an image of the size of the largest input image.
        If the images are asymmetric in opposite directions, will return the
        largest image in both directions.
    return_fft - Return the FFT instead of the convolution.  Useful for making
        PSDs.
    fftshift - If return_fft on, will shift & crop image to appropriate
        dimensions
    kernel_number - if specified, assumes that image # is the kernel, otherwise
        uses the smaller of the two or the first if they are equal-sized
    ignore_nan - attempts to re-weight assuming NAN values are meant to be
        ignored, not treated as zero.  
    ignore_zeros - Ignore the zero-pad-created zeros.  Desirable if you have
        periodic boundaries on a non-2^n grid
    min_wt - If ignoring nans/zeros, force all grid points with a weight less
        than this value to NAN (the weight of a grid point with *no* ignored
        neighbors is 1.0)
    """

    # NAN catching
    nanmask1 = im1!=im1
    nanmask2 = im2!=im2
    im1[nanmask1] = 0
    im2[nanmask2] = 0

    shape1 = im1.shape
    shape2 = im2.shape
    # find ideal size (power of 2) for fft.  Can add shapes because they are tuples
    if fft_pad:
        if psf_pad: 
            # add the X dimensions and Y dimensions and then take the max (bigger)
            fsize = 2**np.ceil(np.log2(np.max(np.array(shape1)+np.array(shape2)))) 
        else: 
            # add the shape lists (max of a list of length 4) (smaller)
            fsize = 2**np.ceil(np.log2(np.max(shape1+shape2)))
        newshape = np.array([fsize,fsize])
    else:
        if psf_pad:
            newshape = np.array(shape1)+np.array(shape2) # just add the biggest dimensions
        else:
            newshape = np.array([np.max([shape1[0],shape2[0]]),np.max([shape1[1],shape2[1]])]) 
    centerx, centery = newshape/2.
    im1quarter1x, im1quarter1y = centerx - shape1[0]/2.,centery - shape1[1]/2.
    im1quarter3x, im1quarter3y = centerx + shape1[0]/2.,centery + shape1[1]/2.
    im2quarter1x, im2quarter1y = centerx - shape2[0]/2.,centery - shape2[1]/2.
    im2quarter3x, im2quarter3y = centerx + shape2[0]/2.,centery + shape2[1]/2.
    bigim1 = np.zeros(newshape,dtype=np.complex128)
    bigim2 = np.zeros(newshape,dtype=np.complex128)
    bigim1[im1quarter1x:im1quarter3x,im1quarter1y:im1quarter3y] = im1
    bigim2[im2quarter1x:im2quarter3x,im2quarter1y:im2quarter3y] = im2 
    imfft1 = fft2(bigim1)
    imfft2 = fft2(bigim2)
    fftmult = imfft1*imfft2
    if ignore_nan or ignore_zeros:
        bigimwt = np.ones(newshape,dtype=np.complex128)
        if ignore_zeros: bigimwt[:] = 0.0
        if shape1[0] > shape2[0] or kernel_number==2: 
            bigimwt[im1quarter1x:im1quarter3x,im1quarter1y:im1quarter3y] = 1.0-nanmask1*ignore_nan
            kern = imfft2
        else:
            bigimwt[im2quarter1x:im2quarter3x,im2quarter1y:im2quarter3y] = 1.0-nanmask2*ignore_nan
            kern = imfft1
        wtfft = fft2(bigimwt)
        wtfftmult = wtfft*kern
        wtfftsm   = ifft2(wtfftmult)
        pixel_weight = np.fft.fftshift(wtfftsm).real

    quarter1x = np.min([im1quarter1x,im2quarter1x])
    quarter3x = np.max([im1quarter3x,im2quarter3x])
    quarter1y = np.min([im1quarter1y,im2quarter1y])
    quarter3y = np.max([im1quarter3y,im2quarter3y])
    if return_fft: 
        if fftshift: # default on 
            if crop:
                return np.fft.fftshift(fftmult)[ quarter1x:quarter3x, quarter1y:quarter3y ]
            else:
                return np.fft.fftshift(fftmult)
        else:
            return fftmult

    if ignore_nan or ignore_zeros:
        rifft = np.fft.fftshift( ifft2( fftmult ) ) / pixel_weight
        rifft[pixel_weight < min_wt] = np.nan
    else:
        rifft = np.fft.fftshift( ifft2( fftmult ) ) 
    if crop:
        result = rifft[ quarter1x:quarter3x, quarter1y:quarter3y ].real
        return result
    else:
        return rifft.real

def smooth(image, kernelwidth=3, kerneltype='gaussian', trapslope=None,
        silent=True, psf_pad=True, interp_nan=False, nwidths='max',
        min_nwidths=6, return_kernel=False, **kwargs):
    """
    Returns a smoothed image using a gaussian, boxcar, or tophat kernel

    Options: 
    kernelwidth - width of kernel in pixels  (see definitions below)
    kerneltype - gaussian, boxcar, or tophat.  
        For a gaussian, uses a gaussian with sigma = kernelwidth (in pixels)
            out to [nwidths]-sigma
        A boxcar is a kernelwidth x kernelwidth square 
        A tophat is a flat circle with radius = kernelwidth
    Default options:
        psf_pad = True - will pad the input image to be the image size + PSF.
            Slows things down but removes edge-wrapping effects (see convolve)
            This option should be set to false if the edges of your image are
            symmetric.
        interp_nan = False - Will replace NaN points in an image with the
            smoothed average of its neighbors
        silent = True - turn it off to get verbose statements about kernel
            types
        return_kernel = False - If set to true, will return the kernel as the
            second return value
        nwidths = 'max' - number of kernel widths wide to make the kernel.  Set to
            'max' to match the image shape, otherwise use any integer 
        min_nwidths = 6 - minimum number of gaussian widths to make the kernel
            (the kernel will be larger than the image if the image size is <
            min_widths*kernelsize)

    Note that the kernel is forced to be even sized on each axis to assure no
    offset when smoothing.
    """

    if (kernelwidth*min_nwidths > image.shape[0] or kernelwidth*min_nwidths > image.shape[1]):
        nwidths = min_nwidths
    if (nwidths!='max'):# and kernelwidth*nwidths < image.shape[0] and kernelwidth*nwidths < image.shape[1]):
        dimsize = np.ceil(kernelwidth*nwidths)
        dimsize += dimsize % 2
        xx,yy = np.indices([dimsize,dimsize])
        sz1,sz2 = dimsize,dimsize
    else:
        sz1,sz2 = image.shape
        sz1 += sz1 % 2
        sz2 += sz2 % 2
        xx,yy = np.indices([sz1,sz2])
    shape = (sz1,sz2)
    if not silent: print "Kernel size set to ",shape
    rr = np.sqrt((xx-sz2/2.)**2+(yy-sz1/2.)**2)

    if kerneltype == 'gaussian':
        kernel = np.exp(-(rr**2)/(2.*kernelwidth**2))
        kernel /= kernel.sum() #/ (kernelwidth**2 * (2*np.pi))
#        if kernelwidth != np.round(kernelwidth):
#            print "Rounding kernel width to %i pixels" % np.round(kernelwidth)
#            kernelwidth = np.round(kernelwidth)

    if kerneltype == 'boxcar':
        if not silent: print "Using boxcar kernel size %i" % np.ceil(kernelwidth)
        kernel = np.ones([np.ceil(kernelwidth),np.ceil(kernelwidth)],dtype='float64') / kernelwidth**2
    elif kerneltype == 'tophat':
        kernel = np.zeros(shape,dtype='float64')
        kernel[rr<kernelwidth] = 1.0
        # normalize
        kernel /= kernel.sum()
    elif kerneltype == 'brickwall':
        if not silent: print "Smoothing with a %i pixel airy function" % kernelwidth
        # airy function is first bessel(x) / x  [like the sinc]
        kernel = j1(rr/kernelwidth) / (rr/kernelwidth) 
        # fix NAN @ center
        kernel[rr==0] = 0.5
        # normalize - technically, this should work, but practically, flux is GAINED in some images.  
        kernel /= kernel.sum()
    elif kerneltype == 'trapezoid':
        if trapslope:
            zz = rr.max()-(rr*trapslope)
            zz[zz<0] = 0
            zz[rr<kernelwidth] = 1.0
            kernel = zz/zz.sum()
        else:
            if not silent: print "trapezoid function requires a slope"

    bad = (image != image)
    temp = image.copy()
    # convolve does this already temp[bad] = 0
    temp = convolve(temp,kernel,kernel_number=2,ignore_nan=interp_nan,psf_pad=psf_pad,**kwargs)
    if interp_nan is False: temp[bad] = image[bad]

    if return_kernel: return temp,kernel
    else: return temp

