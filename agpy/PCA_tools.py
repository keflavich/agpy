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

def PCA_linear_fit(data1, data2):
    """
    Use principal component analysis to determine the best linear fit to the data.
    data1 - x array
    data2 - y array

    returns m,b in the equation y = m x + b
    """
    
    arr = numpy.array([data1-data1.mean(),data2-data2.mean()])

    covmat = numpy.dot(arr,arr.T)
    evals,evects = numpy.linalg.eig(covmat)

    m = evects[1,0] / evects[0,0]
    b = data2.mean() - m*data1.mean()

    return m,b


def smooth_waterfall(arr,fwhm=4.0,unsharp=False):
    """
    Smooth a waterfall plot.

    If unsharp set, remove the smoothed component
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

