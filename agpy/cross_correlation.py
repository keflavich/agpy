"""
Module to measure spectral cross-correlation...
"""
import numpy as np
import lmfit

def fit_lag(arr1,arr2,kind='linear'):

    if arr1.size != arr2.size:
        raise ValueError("Size mismatch")

    def mychi2(pars):
        lag = pars['lag'].value
        return chi2(arr1,arr2,lag)
        
    fitpars = lmfit.Parameters()
    fitpars['lag'] = lmfit.Parameter(value=0.0)

    fitter = lmfit.minimize(mychi2,fitpars,args=())

    return fitter

def shift(data, deltax, phase=0):
    """
    FFT-based sub-pixel image shift
    http://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation/content/html/efficient_subpixel_registration.html
    """
    nx = data.size
    Nx = np.fft.ifftshift(np.linspace(-np.fix(nx/2),np.ceil(nx/2)-1,nx))
    gg = np.fft.ifft( np.fft.fft(data)* np.exp(1j*2*np.pi*(-deltax*Nx/nx)) * np.exp(-1j*phase) )
    return gg

def chi2(arr1,arr2,lag):
    from scipy.interpolate import interp1d
    xv = np.arange(arr1.size)
    #shifted = interp1d(xv, arr2, bounds_error=False, kind='cubic')(xv+lag)
    #shifted[(xv+lag > xv.max()) + (xv+lag < xv.min())] = 0
    #arr1cp = arr1 * (xv>lag) * (xv < (xv+lag).max())
    #shifted = np.interp(xv, xv+lag, arr2, left=np.nan, right=np.nan)
    shifted = np.real(shift(arr2, lag))
    shifted[(xv-lag > xv.max()) + (xv-lag < xv.min())] = 0
    arr1cp = arr1 * (xv>=lag) * (xv <= (xv+lag).max())
    ngood = (shifted != 0).sum()
    if (arr1cp!=0).sum() != ngood:
        print "There are %i OK in input and %i ok in comparison" % ((arr1cp!=0).sum(),ngood)
    if np.any(np.isnan(shifted)):
        raise ValueError("Uncaught NAN")
    return (arr1cp-shifted) / ngood**0.5


if __name__ == "__main__":
    print "Running test code"
    from pylab import *


    xvals = np.linspace(-1,1,100)
    inds = np.arange(xvals.size)
    # gaussian
    test_spec_1 = exp(-xvals**2/(2*0.25**2))
    # power-law
    test_spec_2 = np.abs(np.fft.fftshift( 
        np.fft.fft(np.random.randn(xvals.size)*xvals**-2 +
        np.random.randn(xvals.size)*xvals**-2)*1j ))
    test_spec_2 *= test_spec_1.sum()/test_spec_2.sum() 

    figure(1)
    clf()
    plot(xvals,test_spec_1)
    plot(xvals,test_spec_2)

    # total s/n = 100
    noise = test_spec_1.sum() / 200. 
    #noise = test_spec_1.sum() / 1000. 

    noise1 = np.random.randn(xvals.size)*noise
    test_spec_1n = test_spec_1+noise1
    noise2 = np.random.randn(xvals.size)*noise
    test_spec_2n = test_spec_2+noise2 

    figure(2)
    clf()
    plot(xvals,test_spec_1n)
    plot(xvals,test_spec_2n)

    figure(3)
    clf()
    plot(xvals,np.correlate(test_spec_1n,test_spec_1n,mode='same')/np.correlate(np.ones(xvals.size),np.ones(xvals.size),mode='same'))
    plot(xvals,np.correlate(test_spec_2n,test_spec_2n,mode='same')/np.correlate(np.ones(xvals.size),np.ones(xvals.size),mode='same'))

    figure(4)
    clf()
    plot(xvals,np.correlate(noise1,noise1,mode='same')/np.correlate(np.ones(xvals.size),np.ones(xvals.size),mode='same'))
    plot(xvals,np.correlate(noise2,noise2,mode='same')/np.correlate(np.ones(xvals.size),np.ones(xvals.size),mode='same'))


    chi2_1_wrap = np.array([((test_spec_1-roll(test_spec_1n,ii))**2).sum() for ii in (xvals.size/2-inds)])
    chi2_2_wrap = np.array([((test_spec_2-roll(test_spec_2n,ii))**2).sum() for ii in (xvals.size/2-inds)])

    figure(5)
    clf()
    plot(xvals.size-inds,chi2_1_wrap)
    plot(xvals.size-inds,chi2_2_wrap)
    ylabel("$\\chi^2$")
    xlabel("Lag")

    chi2_1_chop = np.array( 
            [((test_spec_1[:xvals.size-ii]-test_spec_1n[ii:])**2).sum() for ii in xrange(xvals.size/2)]+
            [((test_spec_1[ii:]-test_spec_1n[:xvals.size-ii])**2).sum() for ii in xrange(xvals.size/2,0,-1)]
            )
    chi2_2_chop = np.array( 
            [((test_spec_2[:xvals.size-ii]-test_spec_2n[ii:])**2).sum() for ii in xrange(xvals.size/2)]+
            [((test_spec_2[ii:]-test_spec_2n[:xvals.size-ii])**2).sum() for ii in xrange(xvals.size/2,0,-1)]
            )

    figure(6)
    clf()
    plot(xvals.size/2-inds,chi2_1_chop/(100-np.abs(np.linspace(50,-50,100))))
    plot(xvals.size/2-inds,chi2_2_chop/(100-np.abs(np.linspace(50,-50,100))))
    ylabel("$\\chi^2$")
    xlabel("Lag")


    chi2_smallshifts_1n = [(chi2(test_spec_1,test_spec_1n,xx)**2).sum() for xx in np.linspace(-5,5)]
    chi2_smallshifts_2n = [(chi2(test_spec_2,test_spec_2n,xx)**2).sum() for xx in np.linspace(-5,5)]
    figure(7)
    clf()
    plot(np.linspace(-5,5),chi2_smallshifts_1n)
    plot(np.linspace(-5,5),chi2_smallshifts_2n)
    ylabel("$\\chi^2$")
    xlabel("Lag")
    print "Best-fit shifts (real is zero): ",linspace(-5,5)[argmin(chi2_smallshifts_1n)],linspace(-5,5)[argmin(chi2_smallshifts_2n)]

    from scipy.interpolate import interp1d
    offset_spectra_1n = np.array([interp1d(inds, test_spec_1n, bounds_error=False, kind='cubic')(inds+lag) for lag in linspace(0,1,5)])
    offset_spectra_2n = np.array([interp1d(inds, test_spec_2n, bounds_error=False, kind='cubic')(inds+lag) for lag in linspace(0,1,5)])
    #offset_spectra_1n = np.array([interp(inds, (inds+lag), test_spec_1n) for lag in linspace(0,1,5)])
    #offset_spectra_2n = np.array([interp(inds, (inds+lag), test_spec_2n) for lag in linspace(0,1,5)])
    figure(8)
    clf()
    plot(xvals,offset_spectra_1n.T)
    plot(xvals,offset_spectra_2n.T)

    figure(9)
    clf()
    plot(xvals,(test_spec_1-offset_spectra_1n).T)
    plot(xvals,(test_spec_2-offset_spectra_2n).T)


    import scipy
    import scipy.signal
    xc_1n = scipy.signal.correlate(test_spec_1,test_spec_1n,mode='same')
    xc_2n = scipy.signal.correlate(test_spec_2,test_spec_2n,mode='same')
    figure(10)
    clf()
    plot(xc_1n)
    plot(xc_2n)
    print argmax(xc_1n),argmax(xc_2n)

    shifted_1n = np.interp(inds, inds+1.75, test_spec_1n)
    shifted_2n = np.interp(inds, inds+1.75, test_spec_2n)

    fit1n = fit_lag(test_spec_1,shifted_1n)
    fit2n = fit_lag(test_spec_2,shifted_2n)
    print fit1n.params
    print fit2n.params

    xc_1n = scipy.signal.correlate(test_spec_1,shifted_1n,mode='same')
    xc_2n = scipy.signal.correlate(test_spec_2,shifted_2n,mode='same')
    figure(11)
    clf()
    plot(xc_1n)
    plot(xc_2n)
    print argmax(xc_1n),argmax(xc_2n)

    def fftcorr(s1,s2,pad=0):
        fs1 = np.fft.fft(s1-s1.mean())
        fs2 = np.fft.fft(s2[::-1]-s2.mean())
        mult = fs1*fs2
        if pad:
            mult = np.concatenate([np.zeros(pad),mult,np.zeros(pad)])
        xc = np.fft.ifft(mult)
        return np.fft.fftshift(np.abs(xc))


    test_spec_1_padded = np.concatenate([np.zeros(300), test_spec_1, np.zeros(300)])
    test_spec_2_padded = np.concatenate([np.zeros(300), test_spec_2, np.zeros(300)])
    shifted_1n_padded = np.concatenate([np.zeros(300), shifted_1n, np.zeros(300)])
    shifted_2n_padded = np.concatenate([np.zeros(300), shifted_2n, np.zeros(300)])
    xc_1np = fftcorr(test_spec_1,shifted_1n,pad=300)
    xc_2np = fftcorr(test_spec_2,shifted_2n,pad=300)
    figure(12)
    clf()
    plot(xc_1np)
    plot(xc_2np)
    print argmax(xc_1np),argmax(xc_2np)

    figure(13)
    clf()
    plot(xvals, test_spec_1)
    plot(xvals, shift(test_spec_1,5.5))
    plot(xvals, test_spec_1-shift(test_spec_1,5.5))
    plot(xvals, shift(test_spec_1,55))

    show()
