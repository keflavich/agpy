from AG_fft_tools import correlate2d
from agpy import gaussfitter
import warnings
import numpy as np

def cross_correlation_shifts_FITS(fitsfile1, fitsfile2, return_cropped_images=False, quiet=True, sigma_cut=False, **kwargs):
    """
    Determine the shift between two FITS images using the cross-correlation
    technique.  Requires montage or hcongrid.

    Parameters
    ----------
    fitsfile1: str
        Reference fits file name
    fitsfile2: str
        Offset fits file name
    return_cropped_images: bool
        Returns the images used for the analysis in addition to the measured
        offsets
    quiet: bool
        Silence messages?
    sigma_cut: bool or int
        Perform a sigma-cut before cross-correlating the images to minimize
        noise correlation?
    """
    import montage
    try:
        import astropy.io.fits as pyfits
        import astropy.wcs as pywcs
    except ImportError:
        import pyfits
        import pywcs
    import tempfile

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
    try:
        xoff_wcs,yoff_wcs = np.inner( np.array([[xoff,0],[0,yoff]]), wcs.wcs.cd )[[0,1],[0,1]]
    except AttributeError:
        xoff_wcs,yoff_wcs = 0,0

    if return_cropped_images:
        return xoff,yoff,xoff_wcs,yoff_wcs,image1,image2_projected
    else:
        return xoff,yoff,xoff_wcs,yoff_wcs
    

def cross_correlation_shifts(image1, image2, errim1=None, errim2=None,
        maxoff=None, verbose=False, gaussfit=False, return_error=False,
        zeromean=True, **kwargs):
    """ Use cross-correlation and a 2nd order taylor expansion to measure the
    offset between two images

    Given two images, calculate the amount image2 is offset from image1 to
    sub-pixel accuracy using 2nd order taylor expansion.

    Parameters
    ----------
    image1: np.ndarray
        The reference image
    image2: np.ndarray
        The offset image.  Must have the same shape as image1
    errim1: np.ndarray [optional]
        The pixel-by-pixel error on the reference image
    errim2: np.ndarray [optional]
        The pixel-by-pixel error on the offset image.  
    maxoff: int
        Maximum allowed offset (in pixels).  Useful for low s/n images that you
        know are reasonably well-aligned, but might find incorrect offsets due to 
        edge noise
    zeromean : bool
        Subtract the mean from each image before performing cross-correlation?
    verbose: bool
        Print out extra messages?
    gaussfit : bool
        Use a Gaussian fitter to fit the peak of the cross-correlation?
    return_error: bool
        Return an estimate of the error on the shifts.  WARNING: I still don't
        understand how to make these agree with simulations.
        The analytic estimate comes from
        http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
        At high signal-to-noise, the analytic version overestimates the error
        by a factor of about 1.8, while the gaussian version overestimates
        error by about 1.15.  At low s/n, they both UNDERestimate the error.
        The transition zone occurs at a *total* S/N ~ 1000 (i.e., the total
        signal in the map divided by the standard deviation of the map - 
        it depends on how many pixels have signal)

    **kwargs are passed to correlate2d, which in turn passes them to convolve.
    The available options include image padding for speed and ignoring NaNs.

    References
    ----------
    From http://solarmuri.ssl.berkeley.edu/~welsch/public/software/cross_cor_taylor.pro

    Examples
    --------
    >>> import numpy as np
    >>> im1 = np.zeros([10,10])
    >>> im2 = np.zeros([10,10])
    >>> im1[4,3] = 1
    >>> im2[5,5] = 1
    >>> import AG_image_tools
    >>> yoff,xoff = AG_image_tools.cross_correlation_shifts(im1,im2)
    >>> im1_aligned_to_im2 = np.roll(np.roll(im1,int(yoff),1),int(xoff),0)
    >>> assert (im1_aligned_to_im2-im2).sum() == 0
    

    """

    if not image1.shape == image2.shape:
        raise ValueError("Images must have same shape.")

    if zeromean:
        image1 = image1 - (image1[image1==image1].mean())
        image2 = image2 - (image2[image2==image2].mean())

    quiet = kwargs.pop('quiet') if 'quiet' in kwargs else not verbose
    ccorr = (correlate2d(image1,image2,quiet=quiet,**kwargs) / image1.size)
    # allow for NaNs set by convolve (i.e., ignored pixels)
    ccorr[ccorr!=ccorr] = 0
    if ccorr.shape != image1.shape:
        raise ValueError("Cross-correlation image must have same shape as input images.  This can only be violated if you pass a strange kwarg to correlate2d.")

    ylen,xlen = image1.shape
    xcen = xlen/2-(1-xlen%2) 
    ycen = ylen/2-(1-ylen%2) 

    if ccorr.max() == 0:
        warnings.warn("WARNING: No signal found!  Offset is defaulting to 0,0")
        return 0,0

    if maxoff is not None:
        if verbose: print "Limiting maximum offset to %i" % maxoff
        subccorr = ccorr[ycen-maxoff:ycen+maxoff+1,xcen-maxoff:xcen+maxoff+1]
        ymax,xmax = np.unravel_index(subccorr.argmax(), subccorr.shape)
        xmax = xmax+xcen-maxoff
        ymax = ymax+ycen-maxoff
    else:
        ymax,xmax = np.unravel_index(ccorr.argmax(), ccorr.shape)
        subccorr = ccorr

    if return_error:
        if errim1 is None:
            errim1 = np.ones(ccorr.shape) * image1[image1==image1].std() 
        if errim2 is None:
            errim2 = np.ones(ccorr.shape) * image2[image2==image2].std() 
        eccorr =( (correlate2d(errim1**2, image2**2,quiet=quiet,**kwargs)+
                   correlate2d(errim2**2, image1**2,quiet=quiet,**kwargs))**0.5 
                   / image1.size)
        if maxoff is not None:
            subeccorr = eccorr[ycen-maxoff:ycen+maxoff+1,xcen-maxoff:xcen+maxoff+1]
        else:
            subeccorr = eccorr

    if gaussfit:
        if return_error:
            pars,epars = gaussfitter.gaussfit(subccorr,err=subeccorr,return_all=True)
            exshift = epars[2]
            eyshift = epars[3]
        else:
            pars,epars = gaussfitter.gaussfit(subccorr,return_all=True)
        xshift = maxoff - pars[2] if maxoff is not None else xcen - pars[2]
        yshift = maxoff - pars[3] if maxoff is not None else ycen - pars[3]
        if verbose: 
            print "Gaussian fit pars: ",xshift,yshift,epars[2],epars[3],pars[4],pars[5],epars[4],epars[5]

    else:

        xshift_int = xmax-xcen
        yshift_int = ymax-ycen

        local_values = ccorr[ymax-1:ymax+2,xmax-1:xmax+2]

        d1y,d1x = np.gradient(local_values)
        d2y,d2x,dxy = second_derivative(local_values)

        fx,fy,fxx,fyy,fxy = d1x[1,1],d1y[1,1],d2x[1,1],d2y[1,1],dxy[1,1]

        shiftsubx=(fyy*fx-fy*fxy)/(fxy**2-fxx*fyy)
        shiftsuby=(fxx*fy-fx*fxy)/(fxy**2-fxx*fyy)

        xshift = -(xshift_int+shiftsubx)
        yshift = -(yshift_int+shiftsuby)

        # http://adsabs.harvard.edu/abs/2003MNRAS.342.1291Z
        # Zucker error

        if return_error:
            ccorrn = ccorr / eccorr**2 / ccorr.size #/ (errim1.mean()*errim2.mean()) #/ eccorr**2
            exshift = (np.abs(-1 * ccorrn.size * fxx/ccorrn[ymax,xmax] *
                    (ccorrn[ymax,xmax]**2/(1-ccorrn[ymax,xmax]**2)))**-0.5) 
            eyshift = (np.abs(-1 * ccorrn.size * fyy/ccorrn[ymax,xmax] *
                    (ccorrn[ymax,xmax]**2/(1-ccorrn[ymax,xmax]**2)))**-0.5) 
            if np.isnan(exshift):
                raise ValueError("Error: NAN error!")

    if return_error:
        return xshift,yshift,exshift,eyshift
    else:
        return xshift,yshift

def second_derivative(image):
    """
    Compute the second derivative of an image
    The derivatives are set to zero at the edges

    Parameters
    ----------
    image: np.ndarray

    Returns
    -------
    d/dx^2, d/dy^2, d/dxdy
    All three are np.ndarrays with the same shape as image.

    """
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

def chi2_shift(im1, im2, err=None, mode='wrap', maxoff=None, **kwargs):
    """
    Determine the best fit offset using `scipy.ndimage.map_coordinates` to
    shift the offset image.

    Parameters
    ----------
        im1
        im2
        mode : 'wrap','constant','reflect','nearest'
            Option to pass to map_coordinates for determining what to do with
            shifts outside of the boundaries.  
        maxoff : None or int
            If set, crop the data after shifting before determining chi2
            (this is a good thing to use; not using it can result in weirdness
            involving the boundaries)

    """
    #xc = correlate2d(im1,im2, boundary=boundary)
    #ac1peak = (im1**2).sum()
    #ac2peak = (im2**2).sum()
    #chi2 = ac1peak - 2*xc + ac2peak

    if not im1.shape == im2.shape:
        raise ValueError("Images must have same shape.")

    if np.any(np.isnan(im1)):
        im1 = im1.copy()
        im1[im1!=im1] = 0
    if np.any(np.isnan(im2)):
        im2 = im2.copy()
        if err is not None:
            err[im2!=im2] = np.inf
        im2[im2!=im2] = 0

    im1 = im1-im1.mean()
    im2 = im2-im2.mean()
    yy,xx = np.indices(im1.shape)
    ylen,xlen = im1.shape
    xcen = xlen/2-(1-xlen%2) 
    ycen = ylen/2-(1-ylen%2) 
    import scipy.ndimage,scipy.optimize
    def chi2(p, **kwargs):
        xsh,ysh = p
        shifted_img = scipy.ndimage.map_coordinates(im2, [yy+ysh,xx+xsh], mode=mode, **kwargs)
        if maxoff is not None:
            xslice = slice(xcen-maxoff,xcen+maxoff,None)
            yslice = slice(ycen-maxoff,ycen+maxoff,None)
            # divide by sqrt(number of samples) = sqrt(maxoff**2)
            residuals = np.ravel((im1[yslice,xslice]-shifted_img[yslice,xslice])) / maxoff
        else:
            xslice = slice(None)
            yslice = slice(None)
            residuals = np.abs(np.ravel((im1-shifted_img))) / im1.size**0.5
        if err is None:
            return residuals
        else:
            shifted_err = scipy.ndimage.map_coordinates(err, [yy+ysh,xx+xsh], mode=mode, **kwargs)
            return residuals / shifted_err[slice,slice] 

    bestfit,cov,info,msg,ier = scipy.optimize.leastsq(chi2, [0.5,0.5], full_output=1)
    return bestfit[0],bestfit[1],cov[0,0],cov[1,1]



try:
    import pytest
    import itertools
    from scipy import interpolate

    shifts = [1,1.5,-1.25,8.2,10.1]
    sizes = [99,100,101]
    amps = [5.,10.,50.,100.,500.,1000.]
    gaussfits = (True,False)

    def make_offset_images(xsh,ysh,imsize, width=3.0, amp=1000.0, noiseamp=1.0,
            xcen=50, ycen=50):
        image = np.random.randn(imsize,imsize) * noiseamp
        Y, X = np.indices([imsize, imsize])
        X -= xcen
        Y -= ycen
        new_r = np.sqrt(X*X+Y*Y)
        image += amp*np.exp(-(new_r)**2/(2.*width**2))

        tolerance = 3. * 1./np.sqrt(2*np.pi*width**2*amp/noiseamp)

        new_image = np.random.randn(imsize,imsize)*noiseamp + amp*np.exp(-((X-xsh)**2+(Y-ysh)**2)/(2.*width**2))

        return image, new_image, tolerance

    def make_extended(imsize, powerlaw=2.0):
        yy,xx = np.indices((imsize,imsize))
        cen = imsize/2-(1-imsize%2) 
        yy -= cen
        xx -= cen
        rr = (xx**2+yy**2)**0.5
        
        powermap = (np.random.randn(imsize,imsize) * rr**(-powerlaw)+
            np.random.randn(imsize,imsize) * rr**(-powerlaw) * 1j)
        powermap[powermap!=powermap] = 0

        newmap = np.abs(np.fft.fftshift(np.fft.fft2(powermap)))

        return newmap

    def make_offset_extended(img, xsh, ysh, noise=1.0, mode='wrap'):
        import scipy
        yy,xx = np.indices(img.shape,dtype='float')
        yy-=ysh
        xx-=xsh
        noise = np.random.randn(*img.shape)*noise
        newimage = scipy.ndimage.map_coordinates(img+noise, [yy,xx], mode=mode)

        return newimage



    @pytest.mark.parametrize(('xsh','ysh','imsize','gaussfit'),list(itertools.product(shifts,shifts,sizes,gaussfits)))
    def test_shifts(xsh,ysh,imsize,gaussfit):
        image,new_image,tolerance = make_offset_images(xsh, ysh, imsize)
        if gaussfit:
            xoff,yoff,exoff,eyoff = cross_correlation_shifts(image,new_image)
            print xoff,yoff,np.abs(xoff-xsh),np.abs(yoff-ysh),exoff,eyoff
        else:
            xoff,yoff = cross_correlation_shifts(image,new_image)
            print xoff,yoff,np.abs(xoff-xsh),np.abs(yoff-ysh) 
        assert np.abs(xoff-xsh) < tolerance
        assert np.abs(yoff-ysh) < tolerance

    def do_n_fits(nfits, xsh, ysh, imsize, gaussfit=False, maxoff=None,
            return_error=False, **kwargs):
        """
        Test code

        Parameters
        ----------
        nfits : int
            Number of times to perform fits
        xsh : float
            X shift from input to output image
        ysh : float
            Y shift from input to output image
        imsize : int
            Size of image (square)
        """
        offsets = [
            cross_correlation_shifts( 
                *make_offset_images(xsh, ysh, imsize, **kwargs)[:2],
                gaussfit=gaussfit, maxoff=maxoff, return_error=return_error)
            for ii in xrange(nfits)]

        return offsets

    def do_n_extended_fits(nfits, xsh, ysh, imsize,  gaussfit=False,
            maxoff=None, return_error=False, powerlaw=2.0, noise=1.0,
            unsharp_mask=False, smoothfactor=5, chi2=False, **kwargs):
        image = make_extended(imsize, powerlaw=powerlaw)

        if chi2:
            fitfunc = chi2_shift
        else:
            fitfunc = cross_correlation_shifts

        if unsharp_mask:
            from AG_fft_tools import smooth
            offsets = []
            for ii in xrange(nfits):
                inim = image-smooth(image,smoothfactor)
                offim = make_offset_extended(image, xsh, ysh, noise=noise, **kwargs)
                offim -= smooth(offim,smoothfactor)
                offsets.append( fitfunc( inim, offim, gaussfit=gaussfit,
                        maxoff=maxoff, return_error=return_error) )
        else:
            offsets = [
                fitfunc( 
                    image,
                    make_offset_extended(image, xsh, ysh, noise=noise, **kwargs),
                    gaussfit=gaussfit, maxoff=maxoff, return_error=return_error)
                for ii in xrange(nfits)]

        return offsets


    #@pytest.mark.parametrize(('xsh','ysh','imsize','amp','gaussfit'),
    #        list(itertools.product(shifts,shifts,sizes,amps,gaussfits)))
    def run_tests(xsh, ysh, imsize, amp, gaussfit, nfits=1000, maxoff=20):
        fitted_shifts = np.array(do_n_fits(nfits, xsh, ysh, imsize, amp=amp, maxoff=maxoff))
        errors = fitted_shifts.std(axis=0)
        x,y,ex,ey = cross_correlation_shifts(
                *make_offset_images(xsh, ysh, imsize, amp=amp)[:2],
                gaussfit=gaussfit, maxoff=maxoff, return_error=True,
                errim1=np.ones([imsize,imsize]),
                errim2=np.ones([imsize,imsize]))
        print "StdDev: %10.3g,%10.3g  Measured: %10.3g,%10.3g "+\
                " Difference: %10.3g, %10.3g  Diff/Real: %10.3g,%10.3g" % (
                errors[0],errors[1], ex,ey,errors[0]-ex,errors[1]-ey,
                (errors[0]-ex)/errors[0], (errors[1]-ey)/errors[1])

        return errors[0],errors[1],ex,ey

    def plot_tests(nfits=25,xsh=1.75,ysh=1.75, imsize=64, amp=10., **kwargs):
        x,y,ex,ey = np.array(do_n_fits(nfits, xsh, ysh, imsize, amp,
            maxoff=12., return_error=True, **kwargs)).T
        from pylab import *
        plot([xsh],[ysh],'kd',markersize=20)
        errorbar(x,y,xerr=ex,yerr=ey,linestyle='none')

    def plot_extended_tests(nfits=25,xsh=1.75,ysh=1.75, imsize=64, noise=1.0,
            maxoff=12., **kwargs):
        x,y,ex,ey = np.array(do_n_extended_fits(nfits, xsh, ysh, imsize, 
            maxoff=maxoff, return_error=True, noise=noise, **kwargs)).T
        print x,y
        from pylab import *
        plot([xsh],[ysh],'kd',markersize=20)
        errorbar(x,y,xerr=ex,yerr=ey,linestyle='none')

    def determine_error_offsets():
        """
        Experiment to determine how wrong the error estimates are
        (WHY are they wrong?  Still don't understand)
        """
        # analytic
        A = np.array([run_tests(1.5,1.5,50,a,False,nfits=200) for a in np.logspace(1.5,3,30)]);
        G = np.array([run_tests(1.5,1.5,50,a,True,nfits=200) for a in np.logspace(1.5,3,30)]);
        print "Analytic offset: %g" % (( (A[:,3]/A[:,1]).mean() + (A[:,2]/A[:,0]).mean() )/2. )
        print "Gaussian offset: %g" % (( (G[:,3]/G[:,1]).mean() + (G[:,2]/G[:,0]).mean() )/2. )
        

except ImportError:
    pass
