"""
Test of agpy.gaussfitter.gaussfit
"""
from pylab import *
import numpy
import agpy

# test 1 on figure 1
figure(1)
clf()

# set up grid
xx,yy = indices([100,100])

# make a 2D gaussian with offset & rotation
gg = agpy.gaussfitter.twodgaussian([0,1,40,60,10,20,30])(xx,yy)
subplot(221)
title('2d gaussian')
imshow(gg)

# Make a new gaussiant with some noise added
gg = agpy.gaussfitter.twodgaussian([0,3,40,60,10,20,30])(xx,yy) + randn(100,100)
subplot(222)
title('2d gaussian with noise')
imshow(gg)

# fit the gaussian
myfitpars = agpy.gaussfitter.gaussfit(gg)
# generate the best-fit 2D gaussian
myfit = agpy.gaussfitter.twodgaussian(myfitpars)(xx,yy)
subplot(223)
title('Fit to noisy gaussian')
imshow(myfit)

# should just be noise
subplot(224)
title('noisy 2d gaussian minus best fit')
imshow(gg-myfit)

# prove that it's just gaussian noise
figure(2)
clf()
title('Histogram of noisy gaussian - best fit')
hist((gg-myfit).ravel(),bins=linspace(-4,4,20))

print "Input parameters:  ",[0,3,40,60,10,20,30]
print "Fitted parameters: ",myfitpars

# redo test 1 on figure 4 with masked array
figure(4)
clf()

# set up grid
xx,yy = indices([100,100])

# Mask Indices
n_masked = 1000
mask_xx,mask_yy = floor(random(n_masked)*100).astype('int'),floor(random(n_masked)*100).astype('int')

# make a 2D gaussian with offset & rotation
ggm = numpy.ma.array(agpy.gaussfitter.twodgaussian([0,1,40,60,10,20,30])(xx,yy),mask=zeros([100,100],dtype='bool'))
ggm.mask[mask_yy,mask_xx] = True
subplot(221)
title('2d gaussian')
imshow(ggm)

# Make a new gaussiant with some noise added
ggmn = numpy.ma.array(agpy.gaussfitter.twodgaussian([0,3,40,60,10,20,30])(xx,yy),mask=zeros([100,100],dtype='bool')) + randn(100,100)
ggmn.mask[mask_yy,mask_xx] = True
subplot(222)
title('2d gaussian with noise')
imshow(ggm)

# fit the gaussian
myfitpars = agpy.gaussfitter.gaussfit(ggmn)
# generate the best-fit 2D gaussian
myfit = agpy.gaussfitter.twodgaussian(myfitpars)(xx,yy)
subplot(223)
title('Fit to noisy gaussian')
imshow(myfit)

# should just be noise
subplot(224)
title('noisy 2d gaussian minus best fit')
imshow(ggmn-myfit)

# prove that it's just gaussian noise
figure(2)
clf()
title('Histogram of noisy gaussian - best fit')
hist((gg-myfit).ravel(),bins=linspace(-4,4,20))

print "Input parameters:  ",[0,3,40,60,10,20,30]
print "Fitted parameters: ",myfitpars

#  Test using the scipy FittingData tutorial
#  http://www.scipy.org/Cookbook/FittingData

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

# Create the gaussian data
Xin, Yin = mgrid[0:201, 0:201]
data = gaussian(3, 100, 100, 20, 40)(Xin, Yin) + random(Xin.shape)

figure(3)
clf()
title('scipy FittingData test')
imshow(data, cmap=cm.gist_earth_r)

params = agpy.gaussfitter.gaussfit(data)
fit = agpy.gaussfitter.twodgaussian(params)

contour(fit(*indices(data.shape)), cmap=cm.copper)
ax = gca()
(height, amplitude, x, y, width_x, width_y, angle) = params

text(0.95, 0.05, """
x : %.1f
y : %.1f
width_x : %.1f
width_y : %.1f""" %(x, y, width_x, width_y),
        fontsize=16, horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes)

show()
    
