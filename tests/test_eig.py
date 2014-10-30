import timeit

import numpy as np
import scipy.linalg as scla
import numpy.linalg as npla


if __name__ == "__main__":
    print " ".join(["%16s" % n for n in ("n","sp","np","speigh","npeigh",'spneighn3')])
    methods = ('scla.eig','npla.eig','scla.eigh','npla.eigh','spneighn3')
    alltimes = {k:{} for k in methods}
    for size in (10,100,200,300,400,500,600,700,800,900,1000):

        #array = np.random.random([2**ii]*ndims)
        setup=("import scipy.linalg as scla; import numpy as np; import numpy.linalg as npla;"
               "array = np.random.randn(%i,%i);" % (size,size) + 
               "spneighn3 = lambda x: scla.eigh(x,eigvals=({sz}-3,{sz}-1))".format(sz=size))

        #print array.mean(),array.std()
        #print [ffttype for ffttype in ('sp','np','w')]
        #print min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % 'sp',setup=setup).repeat(3,100))
        for eig in methods:
            alltimes[eig][size] = (min(timeit.Timer(stmt="%s(array)" % eig,setup=setup).repeat(3,10)))

        print "%16i:" % (int(size)) + \
                "".join(["%17f" % alltimes[x][size] for x in methods])

    for eig in methods:

        sizes = alltimes[eig].keys()
        times = alltimes[eig].values()

        m, b = np.polyfit(np.log10(sizes), np.log10(times), 1)

        print "{0} goes as O(n^{1:0.1f}) with offset {2}".format(eig, m, 10**b)


