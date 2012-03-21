import numpy as np
import types
from AG_image_tools import downsample as downsample_2d

try: 
    import fftw3
    has_fftw = True
    def fftwn(array,nthreads=1):
        array = array.astype('complex')
        outarray = array.copy()
        fft_forward = fftw3.Plan(array,outarray, direction='forward', flags=['estimate'],nthreads=nthreads)
        fft_forward()
        return outarray
    def ifftwn(array,nthreads=1):
        array = array.astype('complex')
        outarray = array.copy()
        fft_backward = fftw3.Plan(array,outarray, direction='backward', flags=['estimate'],nthreads=nthreads)
        fft_backward()
        return outarray
except ImportError:
    fftn = np.fft.fftn
    ifftn = np.fft.ifftn
    has_fftw = False
    # I performed some fft speed tests and found that scipy is slower than numpy
    # http://code.google.com/p/agpy/source/browse/trunk/tests/test_ffts.py

def convolve(img, kernel, crop=True, return_fft=False, fftshift=True,
        fft_pad=True, psf_pad=False, ignore_nan=False, quiet=False,
        ignore_zeros=True, min_wt=1e-8, force_ignore_zeros_off=False,
        normalize_kernel=np.sum, debug=False, nthreads=1):
    """
    Convolve an image with a kernel.  Returns a convolved image with shape =
    img.shape.  Assumes image & kernel are centered.

    .. note:: Order matters; the kernel should be second.

    Parameters
    ----------
    - *img* : An n-dimensional `~np.ndarray` object.  Will be convolved
      with *kernel*
    - *kernel* : An n-dimensional `~np.ndarray` object.  Will be normalized
      if *normalize_kernel* is set.  Assumed to be centered (i.e., shifts
      may result if your kernel is asymmetric)

    Options
    -------
    - *fft_pad* :  
        Default on.  Zero-pad image to the nearest 2^n
    - *psf_pad* :  
        Default off.  Zero-pad image to be at least the sum of the image sizes
        (in order to avoid edge-wrapping when smoothing)
    - *crop* :  
        Default on.  Return an image of the size of the largest input image.
        If the images are asymmetric in opposite directions, will return the
        largest image in both directions.
    - *return_fft* : 
        Return the FFT instead of the convolution.  Useful for making PSDs.
    - *fftshift* :
        If return_fft on, will shift & crop image to appropriate dimensions
    - *ignore_nan* :
        attempts to re-weight assuming NAN values are meant to be ignored, not
        treated as zero.  
    - *ignore_zeros* :
        Ignore the zero-pad-created zeros.  Desirable if you have periodic
        boundaries on a non-2^n grid
    - *force_ignore_zeros_off* :
        You can choose to turn off the ignore-zeros when padding, but this is
        only recommended for purely debug purposes
    - *min_wt* :  
        If ignoring nans/zeros, force all grid points with a weight less than
        this value to NAN (the weight of a grid point with *no* ignored
        neighbors is 1.0)
    - *normalize_kernel* : 
        if specified, function to divide kernel by to normalize it
    - *nthreads* :
        if fftw3 is installed, can specify the number of threads to allow FFTs
        to use.  Probably only helpful for large arrays

    Returns
    -------
    *default* : *img* convolved with *kernel*
    if *return_fft* : fft(*img*) * fft(*kernel*)
      - if *fftshift* : Determines whether the fft will be shifted before returning
    if *crop* == False : Returns the image, but with the fft-padded size
      instead of the input size

    """
    
    # replace fftn if has_fftw so that nthreads can be passed
    global fftn,ifftn
    if has_fftw:
        def fftn(*args, **kwargs):
            return fftwn(*args, nthreads=nthreads, **kwargs)

        def ifftn(*args, **kwargs):
            return ifftwn(*args, nthreads=nthreads, **kwargs)


    # mask catching
    if hasattr(img,'mask'):
        mask = img.mask
        img = np.array(img)
        img[mask] = np.nan
    if hasattr(kernel,'mask'):
        mask = kernel.mask
        kernel = np.array(kernel)
        kernel[mask] = np.nan

    # NAN catching
    nanmaskimg = img!=img
    img[nanmaskimg] = 0
    nanmaskkernel = kernel!=kernel
    kernel[nanmaskkernel] = 0
    if (nanmaskimg.sum() > 0 or nanmaskkernel.sum() > 0) and not ignore_nan and not quiet:
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
            newshape = np.array([np.max([imsh,kernsh]) for imsh,kernsh in zip(imgshape,kernshape)]) 

    # separate each dimension by the padding size...
    # this is to determine the appropriate slice size to get back to the input dimensions
    imgslices = []
    kernslices = []
    for ii,(newdimsize,imgdimsize,kerndimsize) in enumerate(zip(newshape,imgshape,kernshape)):
        center = dimsize/2.
        imgslices += [slice(center - imgdimsize/2., center + imgdimsize/2.)]
        kernslices += [slice(center - kerndimsize/2., center + kerndimsize/2.)]

    bigimg = np.zeros(newshape,dtype=np.complex128)
    bigkernel = np.zeros(newshape,dtype=np.complex128)
    bigimg[imgslices] = img
    bigkernel[kernslices] = kernel 
    imgfft = fftn(bigimg)
    kernfft = fftn(bigkernel)
    fftmult = imgfft*kernfft
    if debug: print "Kernel sum: %g  Bigkernel sum: %g" % (kernel.sum(), bigkernel.sum())
    if ignore_nan or ignore_zeros:
        bigimwt = np.ones(newshape,dtype=np.complex128)
        if ignore_zeros: bigimwt[:] = 0.0
        bigimwt[imgslices] = 1.0-nanmaskimg*ignore_nan
        wtfft = fftn(bigimwt)
        wtfftmult = wtfft*kernfft/kernel.sum()
        wtfftsm   = ifftn(wtfftmult)
        pixel_weight = np.fft.fftshift(wtfftsm).real
        if debug: print "Pixel weight: %g  bigimwt: %g wtfft sum: %g kernfft sum: %g wtfftsm sum: %g" % (pixel_weight.sum(), bigimwt.sum(), wtfft.sum(), kernfft.sum(), wtfftsm.sum())

    if np.isnan(fftmult).any():
            raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    if return_fft: 
        if fftshift: # default on 
            if crop:
                return np.fft.fftshift(fftmult)[ imgslices ]
            else:
                return np.fft.fftshift(fftmult)
        else:
            return fftmult

    if ignore_nan or ignore_zeros:
        rifft = np.fft.fftshift( ifftn( fftmult ) ) / pixel_weight
        rifft[pixel_weight < min_wt] = np.nan
    else:
        rifft = np.fft.fftshift( ifftn( fftmult ) ) 
    if crop:
        result = rifft[ imgslices ].real
        return result
    else:
        return rifft.real
