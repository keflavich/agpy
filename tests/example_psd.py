import agpy.psds as psds
import pyfits
import numpy
import pylab

# grab example PSF image
PSFfile = pyfits.open("PSF.fits")
# pixels are 7.2", but to be general, grab them from the header
cd = numpy.abs(PSFfile[0].header['CD1_1']*3600.0)
psf = PSFfile[0].data

xx,yy = numpy.indices(psf.shape)
rr = numpy.sqrt((xx-psf.shape[0]/2.)**2+(yy-psf.shape[1]/2.)**2)

# generate a normalized 2D gaussian for comparison
# best-fit parameters derived using agpy.gaussfitter
gg = numpy.exp(-rr**2/(2.0*(2.20)**2)) / numpy.exp(-rr**2/(2.0*(2.23)**2)).sum()

freqG,psdG = psds.PSD2(gg)
freqP,psdP = psds.PSD2(psf)

# frequencies are stored as normalized frequencies
pylab.figure(1)
pylab.clf()
pylab.semilogy(freqG/cd,psdG,label='Gaussian')
pylab.semilogy(freqP/cd,psdP,label='PSF')
pylab.legend(loc='best')
pylab.axis([0,0.14,1e-7,1])
pylab.xlabel("Spatial Frequency (1/\")")
pylab.ylabel("Normalized Power")
pylab.savefig("example_psd_frequency.png")

pylab.figure(2)
pylab.clf()
pylab.loglog(cd/freqG,psdG,label='Gaussian')
pylab.loglog(cd/freqP,psdP,label='PSF')
pylab.legend(loc='best')
pylab.axis([7,400,1e-7,1])
pylab.xlabel("Spatial Scale (\")")
pylab.ylabel("Normalized Power")
pylab.savefig("example_psd_scale.png")

#pylab.show()
