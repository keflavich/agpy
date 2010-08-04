"""
Test of agpy.gaussfitter.gaussfit
"""
from pylab import *
import agpy

# set up grid
xx,yy = indices([100,100])

# make a 2D gaussian with offset & rotation
gg = agpy.gaussfitter.twodgaussian([0,1,40,60,10,20,30])(xx,yy)
imshow(gg)

# Make a new gaussiant with some noise added
gg = agpy.gaussfitter.twodgaussian([0,3,40,60,10,20,30])(xx,yy) + randn(100,100)
imshow(gg)

# fit the gaussian
myfitpars = agpy.gaussfitter.gaussfit(gg)
# generate the best-fit 2D gaussian
myfit = agpy.gaussfitter.twodgaussian(myfitpars)(xx,yy)
imshow(myfit)

# should just be noise
imshow(gg-myfit)

# prove that it's just gaussian noise
hist((gg-myfit).ravel(),bins=linspace(-4,4,20))

print "Input parameters:  ",[0,3,40,60,10,20,30]
print "Fitted parameters: ",myfitpars
