from agpy import convolve,smooth
from pylab import *

test_image = zeros([512,512])

# point
test_image[204,204] = 1
# box
test_image[280:297,280:297] = 1/256.

figure(0)
clf()
smoothed_gp8_sm25 ,kernel_gp8_25  = smooth(test_image,25.5,'gaussian',nwidths=8,return_kernel=True)
smoothed_gp9_sm25 ,kernel_gp9_25  = smooth(test_image,25.5,'gaussian',nwidths=9,return_kernel=True)
smoothed_gp10_sm25,kernel_gp10_25 = smooth(test_image,25.5,'gaussian',nwidths=10,return_kernel=True)
smoothed_gpmax_sm25,kernel_gpmax_25 = smooth(test_image,25.5,'gaussian',nwidths='max',return_kernel=True)
subplot(221)
imshow(log10(smoothed_gp9_sm25))
colorbar()
subplot(222)
imshow(log10(smoothed_gpmax_sm25))
colorbar()
subplot(223)
title("9-max diff")
imshow(log10(abs(smoothed_gp9_sm25-smoothed_gpmax_sm25)))
colorbar()
subplot(224)
title("9-8 diff")
imshow(log10(abs(smoothed_gp9_sm25-smoothed_gp8_sm25)))
colorbar()
print "Location of the maximum pixel (should be the same for each image) and the shape of the kernel (should be even)"
print "original: ",argmax(test_image)
print "gp8: " ,argmax(smoothed_gp8_sm25), kernel_gp8_25.shape
print "gp9: " ,argmax(smoothed_gp9_sm25), kernel_gp9_25.shape
print "gp10: ",argmax(smoothed_gp10_sm25),kernel_gp10_25.shape
print "gpmax: ",argmax(smoothed_gpmax_sm25),kernel_gpmax_25.shape
print "Maximum value of the kernel (should be ~same) and number of pixels equal to the maximum (should be 1)"
print "npeak gp8: ", (kernel_gp8_25.max()) ,sum((kernel_gp8_25.max()) ==kernel_gp8_25)
print "npeak gp9: ", (kernel_gp9_25.max()) ,sum((kernel_gp9_25.max()) ==kernel_gp9_25)
print "npeak gp10: ",(kernel_gp10_25.max()),sum((kernel_gp10_25.max())==kernel_gp10_25)
print "npeak gpmax: ",(kernel_gpmax_25.max()),sum((kernel_gpmax_25.max())==kernel_gpmax_25)

for ii,smoothsize in enumerate([10,100]): #20,50,100,128]):

    figure(ii+1)
    clf()
    subplot(221)
    imshow(smooth(test_image,smoothsize,'brickwall',silent=False))
    title('Brickwall Filter (Airy pattern)')
    colorbar()
    subplot(222)
    imshow(smooth(test_image,smoothsize,'gaussian',silent=False))
    title('Gaussian')
    colorbar()
    subplot(223)
    imshow(smooth(test_image,smoothsize,'tophat',silent=False))
    title('Tophat')
    colorbar()
    subplot(224)
    imshow(smooth(test_image,smoothsize,'boxcar',silent=False))
    title('Boxcar')
    colorbar()

    print "smoothsize: ",smoothsize
    print "brickwall sum: ", smooth(test_image,smoothsize,'brickwall').sum()
    print "gaussian sum: ", smooth(test_image,smoothsize,'gaussian').sum()
    print "tophat sum: ", smooth(test_image,smoothsize,'tophat').sum()
    print "boxcar sum: ", smooth(test_image,smoothsize,'boxcar').sum()
    print "nofilter sum: ", test_image.sum()

    draw()

