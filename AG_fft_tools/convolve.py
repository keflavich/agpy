import numpy as np
import types
from AG_image_tools import downsample as downsample_2d
from AG_fft_tools import convolvend as convolve

def smooth(image, kernelwidth=3, kerneltype='gaussian', trapslope=None,
        silent=True, psf_pad=True, interp_nan=False, nwidths='max',
        min_nwidths=6, return_kernel=False, normalize_kernel=np.sum,
        downsample=False, downsample_factor=None, **kwargs):
    """
    Returns a smoothed image using a gaussian, boxcar, or tophat kernel

    Options: 
    - *kernelwidth* :
        width of kernel in pixels  (see definitions below)
    - *kerneltype* :
        gaussian, boxcar, or tophat.  
        For a gaussian, uses a gaussian with sigma = kernelwidth (in pixels)
            out to [nwidths]-sigma
        A boxcar is a kernelwidth x kernelwidth square 
        A tophat is a flat circle with radius = kernelwidth
    Default options:
        - *psf_pad* : [True]
            will pad the input image to be the image size + PSF.
            Slows things down but removes edge-wrapping effects (see convolve)
            This option should be set to false if the edges of your image are
            symmetric.
        - *interp_nan* : [False]
            Will replace NaN points in an image with the
            smoothed average of its neighbors (you can still simply ignore NaN 
            values by setting ignore_nan=True but leaving interp_nan=False)
        - *silent* : [True]
            turn it off to get verbose statements about kernel types
        - *return_kernel* : [False]
            If set to true, will return the kernel as the
            second return value
        - *nwidths* : ['max']
            number of kernel widths wide to make the kernel.  Set to 'max' to
            match the image shape, otherwise use any integer 
        - *min_nwidths* : [6]
            minimum number of gaussian widths to make the kernel
            (the kernel will be larger than the image if the image size is <
            min_widths*kernelsize)
        - *normalize_kernel* :
            Should the kernel preserve the map sum (i.e. kernel.sum() = 1)
            or the kernel peak (i.e. kernel.max() = 1) ?  Must be a *function* that can
            operate on a numpy array
        - *downsample* :
            downsample after smoothing?
        - *downsample_factor* :
            if None, default to kernelwidth

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

    if downsample:
        if downsample_factor is None: downsample_factor = kernelwidth
        if return_kernel: return downsample_2d(temp,downsample_factor),downsample_2d(kernel,downsample_factor)
        else: return downsample_2d(temp,downsample_factor)
    else:
        if return_kernel: return temp,kernel
        else: return temp

