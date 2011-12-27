import numpy as np
import types

try: 
    import fftw3
    has_fftw = True
    def fftw2(array,nthreads=1):
        array = array.astype('complex')
        outarray = array.copy()
        fft_forward = fftw3.Plan(array,outarray, direction='forward', flags=['estimate'],nthreads=nthreads)
        fft_forward()
        return outarray
    def ifftw2(array,nthreads=1):
        array = array.astype('complex')
        outarray = array.copy()
        fft_backward = fftw3.Plan(array,outarray, direction='backward', flags=['estimate'],nthreads=nthreads)
        fft_backward()
        return outarray
except ImportError:
    fft2 = np.fft.fft2
    ifft2 = np.fft.ifft2
    has_fftw = False
    # I performed some fft speed tests and found that scipy is slower than numpy
    # http://code.google.com/p/agpy/source/browse/trunk/tests/test_ffts.py
    # try:
    #     #print "Attempting to import scipy.  If you experience a bus error at this step, it is likely because of a bad scipy install"
    #     import scipy
    #     import scipy.fftpack
    #     fft2 = scipy.fftpack.fft2
    #     ifft2 = scipy.fftpack.ifft2
    #     from scipy.special import j1

    # except ImportError:

def convolve(img, kernel, crop=True, return_fft=False, fftshift=True,
        fft_pad=True, psf_pad=False, ignore_nan=False, quiet=False,
        ignore_zeros=True, min_wt=1e-8, force_ignore_zeros_off=False,
        normalize_kernel=np.sum, debug=False, nthreads=1):
    """
    Convolve an image with a kernel.  Returns something the size of an image.
    Assumes image & kernel are centered

    *NOTE* Order matters; the kernel should be second.

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
    ignore_nan - attempts to re-weight assuming NAN values are meant to be
        ignored, not treated as zero.  
    ignore_zeros - Ignore the zero-pad-created zeros.  Desirable if you have
        periodic boundaries on a non-2^n grid
    force_ignore_zeros_off - You can choose to turn off the ignore-zeros when padding,
        but this is only recommended for purely debug purposes
    min_wt - If ignoring nans/zeros, force all grid points with a weight less
        than this value to NAN (the weight of a grid point with *no* ignored
        neighbors is 1.0)
    normalize_kernel - if specified, function to divide kernel by to normalize it

    nthreads - if fftw3 is installed, can specify the number of threads to
        allow FFTs to use.  Probably only helpful for large arrays
    """
    
    # replace fft2 if has_fftw so that nthreads can be passed
    if has_fftw:
        def fft2(*args, **kwargs):
            return fftw2(*args, nthreads=nthreads, **kwargs)

        def ifft2(*args, **kwargs):
            return ifftw2(*args, nthreads=nthreads, **kwargs)


    # mask catching
    if hasattr(img,'mask'):
        mask = img.mask
        img = np.array(img)
        img[mask] = np.nan

    # NAN catching
    nanmaskimg = img!=img
    img[nanmaskimg] = 0
    if nanmaskimg.sum() > 0 and not ignore_nan and not quiet:
        print "Warning: NOT ignoring nan values even though they are present (they are treated as 0)"

    if (psf_pad or fft_pad) and not ignore_zeros and not force_ignore_zeros_off and not quiet:
        print "Warning: when psf_pad or fft_pad are enabled, ignore_zeros is forced on"
        ignore_zeros=True
    elif force_ignore_zeros_off:
        ignore_zeros=False

    if normalize_kernel: # try this.  If a function is not passed, the code will just crash... I think type checking would be better but PEPs say otherwise...
        kernel = kernel / normalize_kernel(kernel)

    if debug: print "Status: ignore_zeros=",ignore_zeros," force_ignore_zeros_off=",force_ignore_zeros_off," psf_pad=",psf_pad," fft_pad=",fft_pad," normalize_kernel=",normalize_kernel

    imgshape = img.shape
    kernshape = kernel.shape
    # find ideal size (power of 2) for fft.  Can add shapes because they are tuples
    if fft_pad:
        if psf_pad: 
            # add the X dimensions and Y dimensions and then take the max (bigger)
            fsize = 2**np.ceil(np.log2(np.max(np.array(imgshape)+np.array(kernshape)))) 
        else: 
            # add the shape lists (max of a list of length 4) (smaller)
            fsize = 2**np.ceil(np.log2(np.max(imgshape+kernshape)))
        newshape = np.array([fsize,fsize])
    else:
        if psf_pad:
            newshape = np.array(imgshape)+np.array(kernshape) # just add the biggest dimensions
        else:
            newshape = np.array([np.max([imgshape[0],kernshape[0]]),np.max([imgshape[1],kernshape[1]])]) 
    centerx, centery = newshape/2.
    imgquarter1x, imgquarter1y = centerx - imgshape[0]/2.,centery - imgshape[1]/2.
    imgquarter3x, imgquarter3y = centerx + imgshape[0]/2.,centery + imgshape[1]/2.
    kernelquarter1x, kernelquarter1y = centerx - kernshape[0]/2.,centery - kernshape[1]/2.
    kernelquarter3x, kernelquarter3y = centerx + kernshape[0]/2.,centery + kernshape[1]/2.
    bigimg = np.zeros(newshape,dtype=np.complex128)
    bigkernel = np.zeros(newshape,dtype=np.complex128)
    bigimg[imgquarter1x:imgquarter3x,imgquarter1y:imgquarter3y] = img
    bigkernel[kernelquarter1x:kernelquarter3x,kernelquarter1y:kernelquarter3y] = kernel 
    imgfft = fft2(bigimg)
    kernfft = fft2(bigkernel)
    fftmult = imgfft*kernfft
    if debug: print "Kernel sum: %g  Bigkernel sum: %g" % (kernel.sum(), bigkernel.sum())
    if ignore_nan or ignore_zeros:
        bigimwt = np.ones(newshape,dtype=np.complex128)
        if ignore_zeros: bigimwt[:] = 0.0
        bigimwt[imgquarter1x:imgquarter3x,imgquarter1y:imgquarter3y] = 1.0-nanmaskimg*ignore_nan
        wtfft = fft2(bigimwt)
        wtfftmult = wtfft*kernfft/kernel.sum()
        wtfftsm   = ifft2(wtfftmult)
        pixel_weight = np.fft.fftshift(wtfftsm).real
        if debug: print "Pixel weight: %g  bigimwt: %g wtfft sum: %g kernfft sum: %g wtfftsm sum: %g" % (pixel_weight.sum(), bigimwt.sum(), wtfft.sum(), kernfft.sum(), wtfftsm.sum())
    if debug: print "img shape:",imgshape," kern shape:",kernshape," newshape:",newshape,\
                    " imgquarter1x,3x,1y,3y:",imgquarter1x,imgquarter3x,imgquarter1y,imgquarter3y, \
                    " kernelquarter1x,3x,1y,3y:",kernelquarter1x,kernelquarter3x,kernelquarter1y,kernelquarter3y  

    if return_fft: 
        if fftshift: # default on 
            if crop:
                return np.fft.fftshift(fftmult)[ imgquarter1x:imgquarter3x, imgquarter1y:imgquarter3y ]
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
        result = rifft[ imgquarter1x:imgquarter3x, imgquarter1y:imgquarter3y ].real
        return result
    else:
        return rifft.real

def smooth(image, kernelwidth=3, kerneltype='gaussian', trapslope=None,
        silent=True, psf_pad=True, interp_nan=False, nwidths='max',
        min_nwidths=6, return_kernel=False, normalize_kernel=np.sum, **kwargs):
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
            smoothed average of its neighbors (you can still simply ignore NaN 
            values by setting ignore_nan=True but leaving interp_nan=False)
        silent = True - turn it off to get verbose statements about kernel
            types
        return_kernel = False - If set to true, will return the kernel as the
            second return value
        nwidths = 'max' - number of kernel widths wide to make the kernel.  Set to
            'max' to match the image shape, otherwise use any integer 
        min_nwidths = 6 - minimum number of gaussian widths to make the kernel
            (the kernel will be larger than the image if the image size is <
            min_widths*kernelsize)
        normalize_kernel - Should the kernel preserve the map sum (i.e. kernel.sum() = 1)
            or the kernel peak (i.e. kernel.max() = 1) ?  Must be a *function* that can
            operate on a numpy array

    Note that the kernel is forced to be even sized on each axis to assure no
    offset when smoothing.
    """

    if (kernelwidth*min_nwidths > image.shape[0] or kernelwidth*min_nwidths > image.shape[1]):
        nwidths = min_nwidths
    if (nwidths!='max'):# and kernelwidth*nwidths < image.shape[0] and kernelwidth*nwidths < image.shape[1]):
        dimsize = np.ceil(kernelwidth*nwidths)
        dimsize += dimsize % 2
        yy,xx = np.indices([dimsize,dimsize])
        szY,szX = dimsize,dimsize
    else:
        szY,szX = image.shape
        szY += szY % 2
        szX += szX % 2
        yy,xx = np.indices([szY,szX])
    shape = (szY,szX)
    if not silent: print "Kernel size set to ",shape
    rr = np.sqrt((xx-szX/2.)**2+(yy-szY/2.)**2)

    if kerneltype == 'gaussian':
        kernel = np.exp(-(rr**2)/(2.*kernelwidth**2))
        kernel /= normalize_kernel(kernel) #/ (kernelwidth**2 * (2*np.pi))
#        if kernelwidth != np.round(kernelwidth):
#            print "Rounding kernel width to %i pixels" % np.round(kernelwidth)
#            kernelwidth = np.round(kernelwidth)

    elif kerneltype == 'boxcar':
        if not silent: print "Using boxcar kernel size %i" % np.ceil(kernelwidth)
        kernel = np.ones([np.ceil(kernelwidth),np.ceil(kernelwidth)],dtype='float64') / kernelwidth**2
    elif kerneltype == 'tophat':
        kernel = np.zeros(shape,dtype='float64')
        kernel[rr<kernelwidth] = 1.0
        # normalize
        kernel /= normalize_kernel(kernel)
    elif kerneltype == 'brickwall':
        if not silent: print "Smoothing with a %i pixel airy function" % kernelwidth
        # airy function is first bessel(x) / x  [like the sinc]
        print "WARNING: I think the Airy should be (2*Besel(x)/x)^2?" # http://en.wikipedia.org/wiki/Airy_disk
        kernel = j1(rr/kernelwidth) / (rr/kernelwidth) 
        # fix NAN @ center
        kernel[rr==0] = 0.5
        # normalize - technically, this should work, but practically, flux is GAINED in some images.  
        kernel /= normalize_kernel(kernel)
    elif kerneltype == 'trapezoid':
        if trapslope:
            zz = rr.max()-(rr*trapslope)
            zz[zz<0] = 0
            zz[rr<kernelwidth] = 1.0
            kernel = zz/zz.sum()
        else:
            if not silent: print "trapezoid function requires a slope"

    if not silent: print "Kernel of type %s normalized with %s has peak %g" % (kerneltype, normalize_kernel, kernel.max())

    bad = (image != image)
    temp = image.copy() # to preserve NaN values
    # convolve does this already temp[bad] = 0

    # kwargs parsing to avoid duplicate keyword passing
    if not kwargs.has_key('ignore_zeros'): kwargs['ignore_zeros']=True
    if not kwargs.has_key('ignore_nan'): kwargs['ignore_nan']=interp_nan

    # No need to normalize - normalization is dealt with in this code
    temp = convolve(temp,kernel,psf_pad=psf_pad, normalize_kernel=False, **kwargs)
    if interp_nan is False: temp[bad] = image[bad]

    if temp.shape != image.shape:
        raise ValueError("Output image changed size; this is completely impossible.")

    if return_kernel: return temp,kernel
    else: return temp

