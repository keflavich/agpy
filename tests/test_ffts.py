import timeit

import scipy.fftpack
spfft2 = scipy.fftpack.fft2
spifft2 = scipy.fftpack.ifft2

import numpy as np
npfft2 = np.fft.fft2
npifft2 = np.fft.ifft2

import fftw3
def wfft2(array,nthreads=4):
    array = array.astype('complex')
    outarray = array.copy()
    fft_forward = fftw3.Plan(array,outarray, direction='forward', flags=['estimate'], nthreads=nthreads)
    fft_forward()
    return outarray
def wifft2(array,nthreads=4):
    array = array.astype('complex')
    outarray = array.copy()
    fft_backward = fftw3.Plan(array,outarray, direction='backward', flags=['estimate'], nthreads=nthreads)
    fft_backward()
    return outarray

if __name__ == "__main__":
    print " ".join(["%17s" % n for n in ("n","sp","np","fftw")])
    for ii in xrange(3,11):
        array = np.random.random([2**ii,2**ii])
        setup="import test_ffts; import numpy as np; array = np.random.random([%i,%i])" % (2**ii,2**ii)

        #print array.mean(),array.std()
        #print [ffttype for ffttype in ('sp','np','w')]
        #print min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % 'sp',setup=setup).repeat(3,100))

        print "%16i:" % (int(2**ii)) + \
                "".join(
                    ["%17f" % (min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % ffttype,setup=setup).repeat(3,100)))
                        for ffttype in ('sp','np','w')]
                )
                
