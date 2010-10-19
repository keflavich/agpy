from agpy import cubes
import numpy as np
import time

parallel_times = []
nonparallel_times = []
sizes = np.array([25,50,100,125,150,200,250,300])

for size in sizes:
    zz,yy,xx = np.indices([size,size,size])
    rr = np.sqrt((size/2)-((zz-(size/2))**2+(yy-(size/2))**2+(xx-(size/2))**2))

    t0 = time.time()

    smcube = cubes.smooth_cube(rr,ignore_nan=True)
    
    t1=time.time()
    print "parallel smooth took %0.2f seconds for size = %i" % (t1-t0,size)
    parallel_times.append(t1-t0)

    smcube = cubes.smooth_cube(rr,ignore_nan=True,parallel=False)

    t2=time.time()
    print "non-parallel smooth took %0.2f seconds for size = %i" % (t2-t1,size)
    nonparallel_times.append(t2-t1)

print parallel_times
print nonparallel_times

proctime = np.zeros(16)
size = 200
zz,yy,xx = np.indices([size,size,size])
rr = np.sqrt((size/2)-((zz-(size/2))**2+(yy-(size/2))**2+(xx-(size/2))**2))

for nprocs in range(1,17):
    t0 = time.time()

    smcube = cubes.smooth_cube(rr,ignore_nan=True,numcores=nprocs)
    
    t1=time.time()
    print "parallel smooth took %0.2f seconds for nproces = %i" % (t1-t0,nprocs)
    proctime[nprocs-1] = t1-t0

import pylab
pylab.figure()
pylab.subplot(211)
pylab.semilogy(sizes,parallel_times   ,linestyle='-',marker='s',label='Parallel')
pylab.semilogy(sizes,nonparallel_times,linestyle='-',marker='o',label='Serial')
pylab.legend(loc='best')
pylab.xlabel('Cube Size')
pylab.ylabel('Execution Time')
pylab.subplot(212)
pylab.plot(sizes,np.array(parallel_times)/np.array(nonparallel_times),linestyle='-',marker='s')
pylab.xlabel('Cube Size')
pylab.ylabel('Parallel/Serial')
pylab.savefig('executiontime_vs_size_cubesmooth.png')

pylab.figure()
pylab.plot(np.arange(1,17),proctime,linestyle='-',marker='s',color='k')
pylab.xlabel('Number of Processors')
pylab.ylabel('Execution Time')
pylab.savefig('executiontime_vs_nprocs.png')

pylab.show()


