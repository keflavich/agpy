from agpy import convolve,smooth
from pylab import *

test_image = zeros([512,512])

# point
test_image[244,244] = 1
# box
test_image[280:284,280:284] = 1

for smoothsize in [10,20,50,100,128]:

    figure()
    clf()
    subplot(221)
    imshow(smooth(test_image,smoothsize,'brickwall'))
    colorbar()
    subplot(222)
    imshow(smooth(test_image,smoothsize,'gaussian'))
    colorbar()
    subplot(223)
    imshow(smooth(test_image,smoothsize,'tophat'))
    colorbar()
    subplot(224)
    imshow(smooth(test_image,smoothsize,'boxcar'))
    colorbar()

    print "smoothsize: ",smoothsize
    print "brickwall sum: ", smooth(test_image,smoothsize,'brickwall').sum()
    print "gaussian sum: ", smooth(test_image,smoothsize,'gaussian').sum()
    print "tophat sum: ", smooth(test_image,smoothsize,'tophat').sum()
    print "boxcar sum: ", smooth(test_image,smoothsize,'boxcar').sum()
    print "nofilter sum: ", test_image.sum()

    show()
