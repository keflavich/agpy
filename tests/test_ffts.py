import timeit

import scipy.fftpack
spfftn = scipy.fftpack.fftn
spifftn = scipy.fftpack.ifftn

import numpy as np
npfftn = np.fft.fftn
npifftn = np.fft.ifftn

import fftw3
def wfftn(array,nthreads=4):
    array = array.astype('complex')
    outarray = array.copy()
    fft_forward = fftw3.Plan(array,outarray, direction='forward', flags=['estimate'], nthreads=nthreads)
    fft_forward()
    return outarray
def wifftn(array,nthreads=4):
    array = array.astype('complex')
    outarray = array.copy()
    fft_backward = fftw3.Plan(array,outarray, direction='backward', flags=['estimate'], nthreads=nthreads)
    fft_backward()
    return outarray

if __name__ == "__main__":
    for ndims in [1,2,3]:
        print "\n%i-dimensional arrays" % ndims
        print " ".join(["%17s" % n for n in ("n","sp","np","fftw")])

        for ii in xrange(3,15-ndims*2):
            #array = np.random.random([2**ii]*ndims)
            setup="import test_ffts; import numpy as np; array = np.random.random([%i]*%i)" % (2**ii,ndims)

            #print array.mean(),array.std()
            #print [ffttype for ffttype in ('sp','np','w')]
            #print min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % 'sp',setup=setup).repeat(3,100))

            print "%16i:" % (int(2**ii)) + \
                    "".join(
                        ["%17f" % (min(timeit.Timer(stmt="test_ffts.%sfftn(array)" % ffttype,setup=setup).repeat(3,10)))
                            for ffttype in ('sp','np','w')]
                    )
                    
