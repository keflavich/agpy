from AG_fft_tools import correlate2d
import numpy as np

def cross_correlation_shifts_FITS(fitsfile1, fitsfile2, return_cropped_images=False, quiet=True, sigma_cut=False, **kwargs):
    """
    Determine the shift between two FITS images using the cross-correlation
    technique.  Requires montage.
    """
    import montage
    import pyfits
    import tempfile
    import pywcs

    header = pyfits.getheader(fitsfile1)
    temp_headerfile = tempfile.NamedTemporaryFile()
    header.toTxtFile(temp_headerfile.name)

    outfile = tempfile.NamedTemporaryFile()
    montage.wrappers.reproject(fitsfile2, outfile.name, temp_headerfile.name, exact_size=True, silent_cleanup=quiet)
    image2_projected = pyfits.getdata(outfile.name)
    image1 = pyfits.getdata(fitsfile1)
    
    outfile.close()
    temp_headerfile.close()

    if image1.shape != image2_projected.shape:
        raise ValueError("montage failed to reproject images to same shape.")

    if sigma_cut:
        corr_image1 = image1*(image1 > image1.std()*sigma_cut)
        corr_image2 = image2_projected*(image2_projected > image2_projected.std()*sigma_cut)
        OK = (corr_image1==corr_image1)*(corr_image2==corr_image2) 
        if (corr_image1[OK]*corr_image2[OK]).sum() == 0:
            print "Could not use sigma_cut of %f because it excluded all valid data" % sigma_cut
            corr_image1 = image1
            corr_image2 = image2_projected
    else:
        corr_image1 = image1
        corr_image2 = image2_projected

    verbose = kwargs.pop('verbose') if 'verbose' in kwargs else not quiet
    xoff,yoff = cross_correlation_shifts(corr_image1, corr_image2, verbose=verbose,**kwargs)
    
    wcs = pywcs.WCS(header)
    xoff_wcs,yoff_wcs = np.inner( np.array([[xoff,0],[0,yoff]]), wcs.wcs.cd )[[0,1],[0,1]]

    if return_cropped_images:
        return xoff,yoff,xoff_wcs,yoff_wcs,image1,image2_projected
    else:
        return xoff,yoff,xoff_wcs,yoff_wcs
    

def cross_correlation_shifts(image1, image2, maxoff=None, verbose=False, **kwargs):
    """
    From http://solarmuri.ssl.berkeley.edu/~welsch/public/software/cross_cor_taylor.pro

    Given two images, calculate the amount image2 is offset from image1 to
    sub-pixel accuracy using 2nd order taylor expansion.

    :maxoff: Maximum allowed offset (in pixels).  Useful for low s/n images that you know
    are reasonably well-aligned, but might find incorrect offsets

    **kwargs are passed to correlate2d, which in turn passes them to convolve.
    The available options include image padding for speed and ignoring NaNs.

    """

    if not image1.shape == image2.shape:
        raise ValueError("Images must have same shape.")

    quiet = kwargs.pop('quiet') if 'quiet' in kwargs else not verbose
    ccorr = correlate2d(image1,image2,quiet=quiet,**kwargs)
    # allow for NaNs set by convolve (i.e., ignored pixels)
    ccorr[ccorr!=ccorr] = 0
    if ccorr.shape != image1.shape:
        raise ValueError("Cross-correlation image must have same shape as input images.  This can only be violated if you pass a strange kwarg to correlate2d.")

    ylen,xlen = image1.shape
    xcen = xlen/2-1 
    ycen = ylen/2-1 

    if maxoff is not None:
        if verbose: print "Limiting maximum offset to %i" % maxoff
        subccorr = ccorr[ycen-maxoff:ycen+maxoff+1,xcen-maxoff:xcen+maxoff+1]
        ymax,xmax = np.nonzero(subccorr == subccorr.max())
        xmax = xmax+xcen-maxoff
        ymax = ymax+ycen-maxoff
    else:
        ymax,xmax = np.nonzero(ccorr == ccorr.max())

    xshift_int = xmax-xcen
    yshift_int = ymax-ycen

    local_values = ccorr[ymax-1:ymax+2,xmax-1:xmax+2]

    d1y,d1x = np.gradient(local_values)
    d2y,d2x,dxy = second_derivative(local_values)

    fx,fy,fxx,fyy,fxy = d1x[1,1],d1y[1,1],d2x[1,1],d2y[1,1],dxy[1,1]

    shiftsubx=(fyy*fx-fy*fxy)/(fxy**2-fxx*fyy)
    shiftsuby=(fxx*fy-fx*fxy)/(fxy**2-fxx*fyy)

    return -(xshift_int+shiftsubx),-(yshift_int+shiftsuby)

def second_derivative(image):
    shift_right = np.roll(image,1,1)
    shift_right[:,0] = 0
    shift_left = np.roll(image,-1,1)
    shift_left[:,-1] = 0
    shift_down = np.roll(image,1,0)
    shift_down[0,:] = 0
    shift_up = np.roll(image,-1,0)
    shift_up[-1,:] = 0

    shift_up_right = np.roll(shift_up,1,1)
    shift_up_right[:,0] = 0
    shift_down_left = np.roll(shift_down,-1,1)
    shift_down_left[:,-1] = 0
    shift_down_right = np.roll(shift_right,1,0)
    shift_down_right[0,:] = 0
    shift_up_left = np.roll(shift_left,-1,0)
    shift_up_left[-1,:] = 0

    dxx = shift_right+shift_left-2*image
    dyy = shift_up   +shift_down-2*image
    dxy=0.25*(shift_up_right+shift_down_left-shift_up_left-shift_down_right)

    return dxx,dyy,dxy

try:
    import pytest
    import itertools
    from scipy import interpolate

    shifts = [1,1.5,-1.25,8.2,10.1]
    sizes = [99,100,101]

    @pytest.mark.parametrize(('xsh','ysh','imsize'),list(itertools.product(shifts,shifts,sizes)))
    def test_shifts(xsh,ysh,imsize):
        width = 3.0
        amp = 1000.0
        noiseamp = 1.0
        image = np.random.randn(imsize,imsize) * noiseamp
        Y, X = np.indices([imsize, imsize])
        X -= 50
        Y -= 50
        new_r = np.sqrt(X*X+Y*Y)
        image += amp*np.exp(-(new_r)**2/(2.*width**2))

        tolerance = 3. * 1./np.sqrt(2*np.pi*width**2*amp/noiseamp)

        new_image = np.random.randn(imsize,imsize)*noiseamp + amp*np.exp(-((X-xsh)**2+(Y-ysh)**2)/(2.*width**2))
        xoff,yoff = cross_correlation_shifts(image,new_image)
        print xoff,yoff,np.abs(xoff-xsh),np.abs(yoff-ysh) 
        assert np.abs(xoff-xsh) < tolerance
        assert np.abs(yoff-ysh) < tolerance

except ImportError:
    pass
