import timeit

import numpy as np
import scipy.linalg as scla
import numpy.linalg as npla


if __name__ == "__main__":
    print " ".join(["%17s" % n for n in ("n","sp","np")])
    for size in (10,20,50,100,200,300,400,500):

        #array = np.random.random([2**ii]*ndims)
        setup=("import scipy.linalg as scla; import numpy as np; import numpy.linalg as npla;"
               "array = np.random.randn(%i,%i)" % (size,size))

        #print array.mean(),array.std()
        #print [ffttype for ffttype in ('sp','np','w')]
        #print min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % 'sp',setup=setup).repeat(3,100))

        print "%16i:" % (int(size)) + \
                "".join(
                    ["%17f" % (min(timeit.Timer(stmt="%s(array)" % eig,setup=setup).repeat(3,10)))
                        for eig in ('scla.eig','npla.eig')]
                )
                
"""
"""
