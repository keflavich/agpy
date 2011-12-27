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
  try:
      arr[arr.mask] = 0
      arr.mask[:] = 0
  except:
      pass
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

def total_least_squares(data1, data2, print_results=False, ignore_nans=True, intercept=True):
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

    """

    if ignore_nans:
        badvals = numpy.isnan(data1) + numpy.isnan(data2)
        if badvals.sum():
            data1 = data1[True-badvals]
            data2 = data2[True-badvals]
    
    if intercept:
        dm1 = data1.mean()
        dm2 = data2.mean()
    else:
        dm1,dm2 = 0,0

    arr = numpy.array([data1-dm1,data2-dm2]).T

    u,s,v = numpy.linalg.svd(arr)

    # v should be sorted.  
    # this solution should be equivalent to v[1,0] / -v[1,1]
    # but I'm using this: http://stackoverflow.com/questions/5879986/pseudo-inverse-of-sparse-matrix-in-python
    m = v[-1,:-1]/-v[-1,-1]

    varfrac = s[0]/s.sum()*100
    if varfrac < 50:
        raise ValueError("ERROR: SVD/TLS Linear Fit accounts for less than half the variance; this is impossible by definition.")

    if intercept:
        b = data2.mean() - m*data1.mean()
        if print_results:
            print "TLS Best fit y = %g x + %g" % (m,b)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*m+b))**2).sum(),data1.shape[0]-2)
        return m,b
    else:
        if print_results:
            print "TLS Best fit y = %g x" % (m)
            print "The fit accounts for %0.3g%% of the variance." % (varfrac)
            print "Chi^2 = %g, N = %i" % (((data2-(data1*m))**2).sum(),data1.shape[0]-1)
        return m


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
    try:
        arr[arr.mask] = 0
        arr.mask[:] = 0
    except:
        pass
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
    try:
        arr[arr.mask] = 0
        arr.mask[:] = 0
    except:
        pass
    covmat = numpy.dot(arr.T,arr)
    evals,evects = numpy.linalg.eig(covmat)
    efuncarr = numpy.dot(arr,evects)
    efuncarr[:,ncomps:] = 0
    return numpy.inner(efuncarr,evects)

if __name__ == "__main__":

    from pylab import *

    xvals = numpy.linspace(0,100,100)
    yvals = numpy.linspace(0,100,100)
    print "First test: m=1, b=0.  Polyfit: ",polyfit(xvals,yvals,1)
    m,b = PCA_linear_fit(xvals,yvals,print_results=True)
    m,b = total_least_squares(xvals,yvals,print_results=True)

    print "Second test: m=-1, b=0. Polyfit: ",polyfit(xvals,yvals*-1,1)
    m,b = PCA_linear_fit(xvals,yvals*-1,print_results=True)
    m,b = total_least_squares(xvals,yvals*-1,print_results=True)

    print "Third test: m=1, b=1. Polyfit: ",polyfit(xvals,yvals+1,1)
    m,b = PCA_linear_fit(xvals,yvals+1,print_results=True)
    m,b = total_least_squares(xvals,yvals+1,print_results=True)

    print "Fourth test: m~1, b~0. Polyfit: ",polyfit(xvals,yvals+random(100),1)
    m,b = PCA_linear_fit(xvals,yvals+random(100),print_results=True)
    m,b = total_least_squares(xvals,yvals+random(100),print_results=True)

    print "Fourth test: m~~1, b~~0. Polyfit: ",polyfit(xvals,yvals+random(100)*50,1)
    m,b = PCA_linear_fit(xvals,yvals+random(100)*50,print_results=True)
    m,b = total_least_squares(xvals,yvals+random(100)*50,print_results=True)

    xr,yr = random(100),random(100)
    print "Fifth test: no linear fit. Polyfit: ",polyfit(xr,yr,1)
    m,b = PCA_linear_fit(xr,yr,print_results=True)
    m,b = total_least_squares(xr,yr,print_results=True)
