"""
=========
PCA Tools
=========

A set of tools for PCA analysis, singular value decomposition,
total least squares, and other linear fitting methods.

Running this code independently tests the fitting functions with different
types of random data.

"""
import numpy

  
def efuncs(arr, return_others=False):
    """
    Determine eigenfunctions of an array for use with
    PCA cleaning
    """
    if hasattr(arr,'filled'):
        arr = arr.filled(0)
    covmat = numpy.dot(arr.T,arr)
    evals,evects = numpy.linalg.eig(covmat)
    efuncarr = numpy.dot(arr,evects)
    if return_others:
        return efuncarr,covmat,evals,evects
    else:
        return efuncarr

def PCA_linear_fit(data1, data2, print_results=False, ignore_nans=True):
    """
    Use principal component analysis to determine the best linear fit to the data.
    data1 - x array
    data2 - y array

    returns m,b in the equation y = m x + b

    print tells you some information about what fraction of the variance is accounted for

    ignore_nans will remove NAN values from BOTH arrays before computing

    Although this works well for the tests below, it fails horrifically on some
    rather well-behaved data sets.  I don't understand why this is, but that's
    why I wrote the total_least_squares SVD code below.
    """

    if ignore_nans:
        badvals = numpy.isnan(data1) + numpy.isnan(data2)
        if badvals.sum():
            data1 = data1[True-badvals]
            data2 = data2[True-badvals]
    
    arr = numpy.array([data1-data1.mean(),data2-data2.mean()])

    covmat = numpy.dot(arr,arr.T)
    evals,evects = numpy.linalg.eig(covmat)

    max_ind = evals.argmax()
    if max_ind:
        evects = evects[::-1,::-1]

    m = evects[1,0] / evects[0,0]
    b = data2.mean() - m*data1.mean()

    varfrac = evals[max_ind]/evals.sum()*100.
    if varfrac < 50:
        raise ValueError("ERROR: PCA Linear Fit accounts for less than half the variance; this is impossible by definition.")

    if print_results:
        print "PCA Best fit y = %g x + %g" % (m,b)
        print "The fit accounts for %0.3g%% of the variance." % (varfrac)
        print "Chi^2 = %g, N = %i" % (((data2-(data1*m+b))**2).sum(),data1.shape[0]-2)

    return m,b

def total_least_squares(data1, data2, data1err=None, data2err=None,
        print_results=False, ignore_nans=True, intercept=True,
        return_error=False, inf=1e10):
    """
    Use Singular Value Decomposition to determine the Total Least Squares linear fit to the data.
    (e.g. http://en.wikipedia.org/wiki/Total_least_squares)
    data1 - x array
    data2 - y array

    if intercept:
        returns m,b in the equation y = m x + b
    else:
        returns m

    print tells you some information about what fraction of the variance is accounted for

    ignore_nans will remove NAN values from BOTH arrays before computing

    Parameters
    ----------
    data1,data2 : np.ndarray
        Vectors of the same length indicating the 'x' and 'y' vectors to fit
    data1err,data2err : np.ndarray or None
        Vectors of the same length as data1,data2 holding the 1-sigma error values

    """

    if ignore_nans:
        badvals = numpy.isnan(data1) + numpy.isnan(data2)
        if data1err is not None:
            badvals += numpy.isnan(data1err)
        if data2err is not None:
            badvals += numpy.isnan(data2err)
        goodvals = True-badvals
        if goodvals.sum() < 2:
            if intercept:
                return 0,0
            else:
                return 0
        if badvals.sum():
            data1 = data1[goodvals]
            data2 = data2[goodvals]

    
    if intercept:
        dm1 = data1.mean()
        dm2 = data2.mean()
    else:
        dm1,dm2 = 0,0

    arr = numpy.array([data1-dm1,data2-dm2]).T

    U,S,V = numpy.linalg.svd(arr, full_matrices=False)

    # v should be sorted.  
    # this solution should be equivalent to v[1,0] / -v[1,1]
    # but I'm using this: http://stackoverflow.com/questions/5879986/pseudo-inverse-of-sparse-matrix-in-python
    M = V[-1,0]/-V[-1,-1]

    varfrac = S[0]/S.sum()*100
    if varfrac < 50:
        raise ValueError("ERROR: SVD/TLS Linear Fit accounts for less than half the variance; this is impossible by definition.")

    # this is performed after so that TLS gives a "guess"
    if data1err is not None or data2err is not None:
        try:
            from scipy.odr import RealData,Model,ODR
        except ImportError:
            raise ImportError("Could not import scipy; cannot run Total Least Squares")

        def linmodel(B,x):
            if intercept:
                return B[0]*x + B[1]
            else:
                return B[0]*x 

        if data1err is not None:
            data1err = data1err[goodvals]
            data1err[data1err<=0] = inf
        if data2err is not None:
            data2err = data2err[goodvals]
            data2err[data2err<=0] = inf

        if any([data1.shape != other.shape for other in (data2,data1err,data2err)]):
            raise ValueError("Data shapes do not match")

        linear = Model(linmodel)
        data = RealData(data1,data2,sx=data1err,sy=data2err)
        B = data2.mean() - M*data1.mean()
        beta0 = [M,B] if intercept else [M]
        myodr = ODR(data,linear,beta0=beta0)
        output = myodr.run()

        if print_results:
            output.pprint()

        if return_error:
            return numpy.concatenate([output.beta,output.sd_beta])
        else:
            return output.beta



    if intercept:
        B = data2.mean() - M*data1.mean()
        if print_results:
            print "TLS Best fit y = %g x + %g" % (M,B)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*M+B))**2).sum(),data1.shape[0]-2)
        return M,B
    else:
        if print_results:
            print "TLS Best fit y = %g x" % (M)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*M))**2).sum(),data1.shape[0]-1)
        return M


def smooth_waterfall(arr,fwhm=4.0,unsharp=False):
    """
    Smooth a waterfall plot.

    If unsharp set, remove the smoothed component

    Input array should have dimensions [timelen, nbolos]
    """

    timelen,nbolos = arr.shape
    kernel = numpy.exp(-numpy.linspace(-timelen/2,timelen/2,timelen)**2/
            (2.0*fwhm/numpy.sqrt(8*numpy.log(2))))
    kernel /= kernel.sum()
    kernelfft = numpy.fft.fft(kernel)
    arrfft = numpy.fft.fft(arr,axis=0)
    arrconv = numpy.fft.fftshift(
            numpy.fft.ifft(arrfft*
            numpy.outer(kernelfft,numpy.ones(nbolos)), 
            axis=0).real,axes=(0,))
    if unsharp:
        return arr-arrconv
    else:
        return arrconv

def pca_subtract(arr,ncomps):
    """
    Compute the eigenfunctions and values of correlated data, then subtract off
    the *ncomps* most correlated components, transform back to the original
    space, and return that.
    """
    if hasattr(arr,'filled'):
        arr = arr.filled(0)
    covmat = numpy.dot(arr.T,arr)
    evals,evects = numpy.linalg.eig(covmat)
    efuncarr = numpy.dot(arr,evects)
    efuncarr[:,0:ncomps] = 0
    return numpy.inner(efuncarr,evects)

def unpca_subtract(arr,ncomps):
    """
    Like pca_subtract, except `keep` the *ncomps* most correlated components
    and reject the others
    """
    if hasattr(arr,'filled'):
        arr = arr.filled(0)
    covmat = numpy.dot(arr.T,arr)
    evals,evects = numpy.linalg.eig(covmat)
    efuncarr = numpy.dot(arr,evects)
    efuncarr[:,ncomps:] = 0
    return numpy.inner(efuncarr,evects)

def pymc_linear_fit(data1, data2, data1err=None, data2err=None,
        print_results=False, intercept=True, nsample=5000, burn=1000,
        thin=10, return_MC=False, guess=None, ignore_nans=True,
        progress_bar=True):
    old_errsettings = numpy.geterr()
    import pymc # pymc breaks np error settings
    numpy.seterr(**old_errsettings)

    if ignore_nans:
        badvals = numpy.isnan(data1) + numpy.isnan(data2)
        if data1err is not None:
            badvals += numpy.isnan(data1err)
        if data2err is not None:
            badvals += numpy.isnan(data2err)
        goodvals = True-badvals
        if goodvals.sum() < 2:
            if intercept:
                return 0,0
            else:
                return 0
        if badvals.sum():
            data1 = data1[goodvals]
            data2 = data2[goodvals]
            if data1err is not None:
                data1err = data1err[goodvals]
            if data2err is not None:
                data2err = data2err[goodvals]


    if guess is None:
        guess = (0,0)

    xmu = pymc.distributions.Uninformative(name='x_observed',value=data1)
    if data1err is None:
        xdata = pymc.distributions.Normal('x', mu=xmu, observed=True,
                value=data1, tau=1, trace=False)
    else:
        xtau = pymc.distributions.Uninformative(name='x_tau',
                value=1.0/data1err**2, observed=True, trace=False)
        xdata = pymc.distributions.Normal('x', mu=xmu, observed=True,
                value=data1, tau=xtau, trace=False)

    d={'slope':pymc.distributions.Uninformative(name='slope', value=guess[0]), }
    if intercept:
        d['intercept'] = pymc.distributions.Uninformative(name='intercept',
                value=guess[1])

        @pymc.deterministic(trace=False)
        def model(x=xdata,slope=d['slope'],intercept=d['intercept']):
            return x*slope+intercept
    else:
        @pymc.deterministic(trace=False)
        def model(x=xdata,slope=d['slope']):
            return x*slope

    d['f'] = model

    if data2err is None:
        ydata = pymc.distributions.Normal('y', mu=model, observed=True,
                value=data2, tau=1, trace=False)
    else:
        ytau = pymc.distributions.Uninformative(name='y_tau',
                value=1.0/data2err**2, observed=True, trace=False)
        ydata = pymc.distributions.Normal('y', mu=model, observed=True,
                value=data2, tau=ytau, trace=False)
    d['y'] = ydata
    
    MC = pymc.MCMC(d)
    MC.sample(nsample,burn=burn,thin=thin,progress_bar=progress_bar)

    MCs = MC.stats()
    m,em = MCs['slope']['mean'],MCs['slope']['standard deviation']
    if intercept: 
        b,eb = MCs['intercept']['mean'],MCs['intercept']['standard deviation']

    if print_results:
        print "MCMC Best fit y = %g x" % (m),
        if intercept: 
            print " + %g" % (b)
        else:
            print ""
        print "m = %g +/- %g" % (m,em)
        if intercept:
            print "b = %g +/- %g" % (b,eb)
        print "Chi^2 = %g, N = %i" % (((data2-(data1*m))**2).sum(),data1.shape[0]-1)

    if return_MC: 
        return MC
    if intercept:
        return m,b
    else:
        return m
        

if __name__ == "__main__":

    from pylab import *

    md,bd = {},{}

    xvals = numpy.linspace(0,100,100)
    yvals = numpy.linspace(0,100,100)
    md['ideal'],bd['ideal'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0},)
    md['ideal']['poly'],bd['ideal']['poly'] = polyfit(xvals,yvals,1)
    md['ideal']['PCA'],bd['ideal']['PCA']   = PCA_linear_fit(xvals,yvals,print_results=True)
    md['ideal']['TLS'],bd['ideal']['TLS']   = total_least_squares(xvals,yvals,print_results=True)
    md['ideal']['pymc'],bd['ideal']['pymc'] = pymc_linear_fit(xvals,yvals,print_results=True)

    md['neg'],bd['neg'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['neg']['poly'],bd['neg']['poly'] = polyfit(xvals,yvals*-1,1)
    md['neg']['PCA'],bd['neg']['PCA']   = PCA_linear_fit(xvals,yvals*-1,print_results=True)
    md['neg']['TLS'],bd['neg']['TLS']   = total_least_squares(xvals,yvals*-1,print_results=True)
    md['neg']['pymc'],bd['neg']['pymc'] = pymc_linear_fit(xvals,yvals*-1,print_results=True)

    md['intercept'],bd['intercept'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['intercept']['poly'],bd['intercept']['poly'] = polyfit(xvals,yvals+1,1)
    md['intercept']['PCA'],bd['intercept']['PCA']   = PCA_linear_fit(xvals,yvals+1,print_results=True)
    md['intercept']['TLS'],bd['intercept']['TLS']   = total_least_squares(xvals,yvals+1,print_results=True)
    md['intercept']['pymc'],bd['intercept']['pymc'] = pymc_linear_fit(xvals,yvals+1,print_results=True)

    md['noise'],bd['noise'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['noise']['poly'],bd['noise']['poly'] = polyfit(xvals,yvals+random(100),1)
    md['noise']['PCA'],bd['noise']['PCA']   = PCA_linear_fit(xvals,yvals+random(100),print_results=True)
    md['noise']['TLS'],bd['noise']['TLS']   = total_least_squares(xvals,yvals+random(100),print_results=True)
    md['noise']['pymc'],bd['noise']['pymc'] = pymc_linear_fit(xvals,yvals+random(100),print_results=True)

    md['highnoise'],bd['highnoise'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['highnoise']['poly'],bd['highnoise']['poly'] = polyfit(xvals,yvals+random(100)*50,1)
    md['highnoise']['PCA'],bd['highnoise']['PCA']   = PCA_linear_fit(xvals,yvals+random(100)*50,print_results=True)
    md['highnoise']['TLS'],bd['highnoise']['TLS']   = total_least_squares(xvals,yvals+random(100)*50,print_results=True)
    md['highnoise']['pymc'],bd['highnoise']['pymc'] = pymc_linear_fit(xvals,yvals+random(100)*50,print_results=True)

    md['random'],bd['random'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    xr,yr = random(100),random(100)
    md['random']['poly'],bd['random']['poly'] = polyfit(xr,yr,1)
    md['random']['PCA'],bd['random']['PCA']   = PCA_linear_fit(xr,yr,print_results=True)
    md['random']['TLS'],bd['random']['TLS']   = total_least_squares(xr,yr,print_results=True)
    md['random']['pymc'],bd['random']['pymc'] = pymc_linear_fit(xr,yr,print_results=True)

    md['xnoise'],bd['xnoise'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['xnoise']['poly'],bd['xnoise']['poly'] = polyfit(xvals+random(100)*5,yvals+random(100)*5,1)
    md['xnoise']['PCA'],bd['xnoise']['PCA']   = PCA_linear_fit(xvals+random(100)*5,yvals+random(100)*5,print_results=True)
    md['xnoise']['TLS'],bd['xnoise']['TLS']   = total_least_squares(xvals+random(100)*5,yvals+random(100)*5,print_results=True)
    md['xnoise']['pymc'],bd['xnoise']['pymc'] = pymc_linear_fit(xvals+random(100)*5,yvals+random(100)*5,print_results=True)

    md['xhighnoise'],bd['xhighnoise'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['xhighnoise']['poly'],bd['xhighnoise']['poly'] = polyfit(xvals+random(100)*50,yvals+random(100)*50,1)
    md['xhighnoise']['PCA'],bd['xhighnoise']['PCA']   = PCA_linear_fit(xvals+random(100)*50,yvals+random(100)*50,print_results=True)
    md['xhighnoise']['TLS'],bd['xhighnoise']['TLS']   = total_least_squares(xvals+random(100)*50,yvals+random(100)*50,print_results=True)
    md['xhighnoise']['pymc'],bd['xhighnoise']['pymc'] = pymc_linear_fit(xvals+random(100)*50,yvals+random(100)*50,print_results=True)

    md['witherrors'],bd['witherrors'] = ({'PCA':0,'TLS':0, 'poly':0, 'pymc':0},{'PCA':0,'TLS':0, 'poly':0, 'pymc':0})
    md['witherrors']['poly'],bd['witherrors']['poly'] = 0,0
    md['witherrors']['PCA'],bd['witherrors']['PCA']   = 0,0
    xerr,yerr = randn(100)/10.0+1.0+sqrt(xvals),randn(100)/10.0+1.0+sqrt(yvals)
    x,y = xvals+xerr,yvals+yerr
    md['witherrors']['TLS'],bd['witherrors']['TLS']   = total_least_squares(x,y,data1err=xerr,data2err=yerr,print_results=True)
    md['witherrors']['pymc'],bd['witherrors']['pymc'] = pymc_linear_fit(x,y,data1err=xerr,data2err=yerr,print_results=True)

    print "Slopes: "
    toprow = (" "*20) + " ".join(["%20s" % k for k in md['ideal']])
    print toprow
    for colname,column in md.items():
        print "%20s" % colname,
        for rowname,row in column.items():
            print "%20s" % row,
        print
            
    print "Intercepts: "
    toprow = (" "*20) + " ".join(["%20s" % k for k in bd['ideal']])
    print toprow
    for colname,column in bd.items():
        print "%20s" % colname,
        for rowname,row in column.items():
            print "%20s" % row,
        print
            

    print "PyMC linear tests"
    MC1 = pymc_linear_fit(x,y,intercept=False,print_results=True,return_MC=True)
    MC2 = pymc_linear_fit(x,y,xerr,yerr,intercept=False,print_results=True,return_MC=True)

    hoggdata = np.array([
        [1,201,592,61,9,-0.84],
        [2,244,401,25,4,0.31],
        [3,47,583,38,11,0.64],
        [4,287,402,15,7,-0.27],
        [5,203,495,21,5,-0.33],
        [6,58,173,15,9,0.67],
        [7,210,479,27,4,-0.02],
        [8,202,504,14,4,-0.05],
        [9,198,510,30,11,-0.84],
        [10,158,416,16,7,-0.69],
        [11,165,393,14,5,0.30],
        [12,201,442,25,5,-0.46],
        [13,157,317,52,5,-0.03],
        [14,131,311,16,6,0.50],
        [15,166,400,34,6,0.73],
        [16,160,337,31,5,-0.52],
        [17,186,423,42,9,0.90],
        [18,125,334,26,8,0.40],
        [19,218,533,16,6,-0.78],
        [20,146,344,22,5,-0.56],
        ])

    linear_fitters = [total_least_squares, PCA_linear_fit, pymc_linear_fit]

    for method in linear_fitters:
        print method.__name__,method(hoggdata[:,1],hoggdata[:,2])
        try:
            print method.__name__,method(hoggdata[:,1],hoggdata[:,2], data2err=hoggdata[:,3])
        except TypeError:
            pass
        except AttributeError:
            pass
        try:
            print method.__name__,method(hoggdata[:,1],hoggdata[:,2], data2err=hoggdata[:,3], data1err=hoggdata[:,4])
        except TypeError:
            pass
