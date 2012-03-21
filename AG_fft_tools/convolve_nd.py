import numpy as np
import types
from AG_image_tools import downsample as downsample_2d

try: 
    import fftw3
    has_fftw = True
    def fftwn(array,nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_forward = fftw3.Plan(array,outarray, direction='forward', flags=['estimate'],nthreads=nthreads)
        fft_forward.execute()
        return outarray
    def ifftwn(array,nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_backward = fftw3.Plan(array,outarray, direction='backward', flags=['estimate'],nthreads=nthreads)
        fft_backward.execute()
        return outarray / np.size(array)
except ImportError:
    fftn = np.fft.fftn
    ifftn = np.fft.ifftn
    has_fftw = False
    # I performed some fft speed tests and found that scipy is slower than numpy
    # http://code.google.com/p/agpy/source/browse/trunk/tests/test_ffts.py

def convolvend(img, kernel, crop=True, return_fft=False, fftshift=True,
        fft_pad=True, psf_pad=False, ignore_nan=False, quiet=False,
        ignore_zeros=True, min_wt=1e-8, force_ignore_zeros_off=False,
        normalize_kernel=np.sum, debug=False, use_numpy_fft=not has_fftw,
        nthreads=1):
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
    if has_fftw and not use_numpy_fft:
        if debug: print "Using FFTW"
        def fftn(*args, **kwargs):
            return fftwn(*args, nthreads=nthreads, **kwargs)

        def ifftn(*args, **kwargs):
            return ifftwn(*args, nthreads=nthreads, **kwargs)
    elif use_numpy_fft:
        if debug: print "Using numpy"
        fftn = np.fft.fftn
        ifftn = np.fft.ifftn


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
    if debug:
        print "nanmasked in image: %i kernel: %i" % (nanmaskimg.sum(), nanmaskkernel.sum())

    if (psf_pad or fft_pad) and not ignore_zeros and not force_ignore_zeros_off and not quiet:
        print "Warning: when psf_pad or fft_pad are enabled, ignore_zeros is forced on"
        ignore_zeros=True
    elif force_ignore_zeros_off:
        ignore_zeros=False

    if normalize_kernel: # try this.  If a function is not passed, the code will just crash... I think type checking would be better but PEPs say otherwise...
        kernel = kernel / normalize_kernel(kernel)

    if debug: print "Status: ignore_zeros=",ignore_zeros," force_ignore_zeros_off=",force_ignore_zeros_off," psf_pad=",psf_pad," fft_pad=",fft_pad," normalize_kernel=",normalize_kernel," return_fft=",return_fft

    imgshape = img.shape
    kernshape = kernel.shape
    ndim = len(img.shape)
    if ndim != len(kernshape):
        raise ValueError("Image and kernel must have same number of dimensions")
    # find ideal size (power of 2) for fft.  Can add shapes because they are tuples
    if fft_pad:
        if psf_pad: 
            # add the X dimensions and Y dimensions and then take the max (bigger)
            fsize = 2**np.ceil(np.log2(np.max(np.array(imgshape)+np.array(kernshape)))) 
        else: 
            # add the shape lists (max of a list of length 4) (smaller)
            fsize = 2**np.ceil(np.log2(np.max(imgshape+kernshape)))
        newshape = np.array([fsize for ii in range(ndim)])
    else:
        if psf_pad:
            newshape = np.array(imgshape)+np.array(kernshape) # just add the biggest dimensions
        else:
            newshape = np.array([np.max([imsh,kernsh]) for imsh,kernsh in zip(imgshape,kernshape)]) 

    if debug: 
        print "newshape: ",newshape," imgshape: ",imgshape," kernshape: ",kernshape

    # separate each dimension by the padding size...
    # this is to determine the appropriate slice size to get back to the input dimensions
    imgslices = []
    kernslices = []
    for ii,(newdimsize,imgdimsize,kerndimsize) in enumerate(zip(newshape,imgshape,kernshape)):
        center = newdimsize/2.
        imgslices += [slice(center - imgdimsize/2., center + imgdimsize/2.)]
        kernslices += [slice(center - kerndimsize/2., center + kerndimsize/2.)]

    if debug: 
        print "imgslices: ",imgslices," kernel slices: ",kernslices

    bigimg = np.zeros(newshape,dtype=np.complex128)
    bigkernel = np.zeros(newshape,dtype=np.complex128)
    bigimg[imgslices] = img
    bigkernel[kernslices] = kernel 
    imgfft = fftn(bigimg)
    kernfft = fftn(bigkernel)
    fftmult = imgfft*kernfft
    if debug: 
        print "Kernel sum: %g  Bigkernel sum: %g image sum: %g bigimg sum: %g" % (kernel.sum(), bigkernel.sum(), img.sum(), bigimg.sum())
        print "fft(img) sum: %g  fft(kernel) sum: %g fft(img)*fft(kernel) sum: %g" % (imgfft.sum(), kernfft.sum(), fftmult.sum())
    if ignore_nan or ignore_zeros:
        if ignore_zeros: 
            bigimwt = np.zeros(newshape,dtype=np.complex128)
        else:
            bigimwt = np.ones(newshape,dtype=np.complex128)
        bigimwt[imgslices] = 1.0-nanmaskimg*ignore_nan
        if debug: print "bigimwt.sum: %g" % (bigimwt.sum())
        wtfft = fftn(bigimwt)
        wtfftmult = wtfft*kernfft/kernel.sum() # I think this one HAS to be normalized
        wtsm   = ifftn(wtfftmult)
        # need to re-zero weights outside of the image (if it is padded, we still don't weight those regions)
        bigimwt[imgslices] = np.fft.fftshift(wtsm).real[imgslices]
        if debug:
            print "Pixel weight: %g  bigimwt: %g wtfft sum: %g kernfft sum: %g wtsm sum: %g " % (bigimwt.sum(), bigimwt.sum(), wtfft.sum(), kernfft.sum(), wtsm.sum())
            print "Nonzero weights: bigimwt %i, wtsm %i" % ((bigimwt>0).sum(), (wtsm>0).sum())

    if np.isnan(fftmult).any():
            raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    # restore nans in original image 
    img[nanmaskimg] = np.nan
    kernel[nanmaskkernel] = np.nan

    if return_fft: 
        if fftshift: # default on 
            if crop:
                return np.fft.fftshift(fftmult)[ imgslices ]
            else:
                return np.fft.fftshift(fftmult)
        else:
            return fftmult

    if ignore_nan or ignore_zeros:
        rifft = np.fft.fftshift( ifftn( fftmult ) ) / bigimwt
        if debug:
            print "%i weights < min_wt = %g" % ((bigimwt < min_wt).sum(), min_wt)
            print "%i VALID weights < min_wt = %g" % ((bigimwt[imgslices] < min_wt).sum(), min_wt)
            print "total weights > min_wt: %g" % (bigimwt[bigimwt>=min_wt].sum())
            print "total VALID weights > min_wt: %g" % (bigimwt[imgslices][bigimwt[imgslices]>=min_wt].sum())
        rifft[bigimwt < min_wt] = np.nan
    else:
        if debug: print ifftn, fftmult.sum()
        rifft = np.fft.fftshift( ifftn( fftmult ) ) 

    if debug:
        print "rifft sum: %g" % (rifft.sum())

    if crop:
        result = rifft[ imgslices ].real
        if debug: print "result sum: %g" % (result.sum())
        return result
    else:
        return rifft.real


import pytest
def test_3d(debug=False):
    img = np.zeros([32,32,32])
    img[15,15,15]=1
    img[15,0,15]=1
    kern = np.zeros([32,32,32])
    kern[14:19,14:19,14:19] = 1

    for psf_pad in (True,False):
        conv1 = convolvend(img, kern, psf_pad=psf_pad, debug=debug)
        conv2 = convolvend(img, kern, psf_pad=psf_pad, use_numpy_fft=True, debug=debug)
        conv4 = convolvend(img, kern, psf_pad=psf_pad, use_numpy_fft=True, force_ignore_zeros_off=True, debug=debug)
        conv3 = convolvend(img, kern, psf_pad=psf_pad, force_ignore_zeros_off=True, debug=debug)

        print "psf_pad=%s" % psf_pad
        print "FFTW, ignore_zeros: %g,%g" % (conv1[15,0,15],conv1[15,15,15])
        print "numpy, ignore_zeros: %g,%g" % (conv2[15,0,15],conv2[15,15,15])
        print "FFTW, no ignore_zeros: %g,%g" % (conv3[15,0,15],conv3[15,15,15])
        print "numpy, no ignore_zeros: %g,%g" % (conv4[15,0,15],conv4[15,15,15])

        
