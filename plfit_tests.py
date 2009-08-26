import cplfit
import plfit
import time
from pylab import *
#from numpy.random import rand
#from numpy import unique,sort,array,asarray,log,sum,min,max,argmin,argmax,arange
import sys
""" 
Test code for fixed, varying xmin.  Alpha is set to 2.5, and argv[1] random power-law
distributions are then fit.  

Timing tests implemented in speedcompare_plfit.py
X=plfit.plexp_inv(rand(ne),1,2.5)
t1=time.time(); p3=plfit.plfit(X,usefortran=False,usecy=True); print time.time()-t1
t1=time.time(); p1=plfit.plfit(X); print time.time()-t1
t1=time.time(); p3=plfit.plfit(X,usefortran=False); print time.time()-t1
"""

if len(sys.argv) > 1:
    ntests = int(sys.argv[1])
    if len(sys.argv) > 2:
        nel = int(sys.argv[2])
        if len(sys.argv) > 3:
            xmin = float(sys.argv[3])
        else:
            xmin = 0.5
    else: 
        nel = 1000
else:
    ntests = 10000

a=zeros(ntests)
for i in range(ntests):
    X=plfit.plexp_inv(rand(nel),xmin,2.5)
    p=plfit.plfit(X,xmin=xmin,quiet=True,silent=True)
    a[i]=p._alpha

h,b = hist(a,bins=30)[:2]
bx = (b[1:]+b[:-1])/2.0

import gaussfitter
p,m,pe,chi2 = gaussfitter.onedgaussfit(bx,h,params=[0,ntests/10.0,2.5,0.5],fixed=[1,0,0,0])

plot(bx,m)

print "XMIN fixed: Alpha = 2.5 (real), %0.3f +/- %0.3f (measured)" % (p[2],p[3])


a=zeros(ntests)
xm=zeros(ntests)
for i in range(ntests):
    X=plfit.plexp_inv(rand(nel),xmin,2.5)
    p=plfit.plfit(X,quiet=True,silent=True)
    a[i]=p._alpha
    xm[i]=p._xmin

figure()
h,b = hist(a,bins=30)[:2]
xlabel('alpha')
bx = (b[1:]+b[:-1])/2.0

p,m,pe,chi2 = gaussfitter.onedgaussfit(bx,h,params=[0,ntests/10.0,2.5,0.5],fixed=[1,0,0,0])
plot(bx,m)
print "XMIN varies: Alpha = 2.5 (real), %0.3f +/- %0.3f (measured)" % (p[2],p[3])

figure()
h,b = hist(xm,bins=30)[:2]
xlabel('xmin')
bx = (b[1:]+b[:-1])/2.0

p,m,pe,chi2 = gaussfitter.onedgaussfit(bx,h,params=[0,ntests/10.0,xmin,0.2],fixed=[1,0,0,0])
plot(bx,m)
print "XMIN varies: XMIN = %0.3f (real), %0.3f +/- %0.3f (measured)" % (xmin,p[2],p[3])
