from agpy import convolve,smooth
from pylab import *

test_image = zeros([512,512])

# point
test_image[204,204] = 1
# box
test_image[280:297,280:297] = 1/256.

shiftest_image = zeros([300,320])
shiftest_image[100:120,180:220]=1.0
print "Testing for a shift"
figure(11)
clf()
smoothed,kernel = smooth(shiftest_image,return_kernel=True)
subplot(221)
title("shiftest_image")
imshow(shiftest_image)
subplot(222)
title("smoothed")
imshow(smoothed)
subplot(223)
title("shiftest_image-smoothed")
imshow(shiftest_image-smoothed)
subplot(224)
title("kernel")
imshow(kernel)


"""
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
"""

"""
print "\n\nDemonstration that you need to ignore_zeros when padding (figure 10)"
figure(10)
clf()
testimage = ones([10,20]) # make a flat image
testimage[5,5] += 1 # add some contrast
smtestimage = smooth(testimage)
smtestimage_nopsfpad = smooth(testimage,psf_pad=False,force_ignore_zeros_off=True)
smtestimage_nofftpad = smooth(testimage,fft_pad=False,force_ignore_zeros_off=True)
smtestimage_ignorenan = smooth(testimage,interp_nan=True,force_ignore_zeros_off=True)
smtestimage_ignorezeros = smooth(testimage,ignore_zeros=True,force_ignore_zeros_off=True)
smtestimage_noz_nopsfpad = smooth(testimage,psf_pad=False,ignore_zeros=True)
smtestimage_noz_nofftpad = smooth(testimage,fft_pad=False,ignore_zeros=True)
smtestimage_noz_nopad = smooth(testimage,fft_pad=False,psf_pad=False,ignore_zeros=True)
smtestimage_nopad = smooth(testimage,fft_pad=False,psf_pad=False,force_ignore_zeros_off=True)
subplot(331)
title("smtestimage_nopad")
imshow(smtestimage_nopad)
subplot(332)
title("smtestimage (default)")
imshow(smtestimage)
subplot(333)
title("smtestimage_nopsfpad")
imshow(smtestimage_nopsfpad)
subplot(334)
title("smtestimage_nofftpad")
imshow(smtestimage_nofftpad)
subplot(335)
title("smtestimage_ignorenan")
imshow(smtestimage_ignorenan)
subplot(336)
title("smtestimage_ignorezeros")
imshow(smtestimage_ignorezeros)
subplot(337)
title("smtestimage_noz_nopad")
imshow(smtestimage_noz_nopad)
subplot(338)
title("smtestimage_noz_nopsfpad")
imshow(smtestimage_noz_nopsfpad)
subplot(339)
title("smtestimage_noz_nofftpad")
imshow(smtestimage_noz_nofftpad)
"""

"""
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

"""
