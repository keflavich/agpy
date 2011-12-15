import agpy
from agpy import psf_fitter,asinh_norm
from pylab import *

if __name__ == "__main__":

    # numerically determine Airy FWHM
    # (see http://en.wikipedia.org/wiki/Airy_disk)
    x = linspace(1.5,1.7,100000)
    airy_fwhm = x[argmin(abs(agpy.psf_fitter._airy_func(x) - 0.5))] * 2

    # numerically determine Gaussian fwhm
    x = linspace(1.1,1.3,100000)
    gaussian_fwhm = x[argmin(abs(agpy.psf_fitter._gaussian_func(x) - 0.5))] * 2
    gaussian_fwhm_analytic = sqrt(8*log(2))

    print "Airy FWHM:     %15.10f" % airy_fwhm
    print "Gaussian FWHM: %15.10f" % gaussian_fwhm
    print "Gaussian FWHM: %15.10f (analytic)" % gaussian_fwhm_analytic
    print "Airy 1-'sigma':%15.10f" % (agpy.psf_fitter._airy_func(1.0))
    print "Gaussian 1sig: %15.10f" % (agpy.psf_fitter._gaussian_func(1.0))

    figure(1)
    clf()
    x = linspace(0,5,1000)
    plot(x,agpy.psf_fitter._airy_func(x),label='Airy')
    plot(x,agpy.psf_fitter._gaussian_func(x),label='Gaussian')
    plot(x,agpy.psf_fitter._gaussian_func(x,sigma=airy_fwhm/gaussian_fwhm),label='Gaussian (sigma=%5.3f)' % (airy_fwhm/gaussian_fwhm))


    figure(2)
    clf()
    test_airy = psf_fitter.airy([1,50,50,5],vheight=False,shape=[100,100])
    imshow((test_airy),norm=asinh_norm.AsinhNorm())
    colorbar()

    airyfitp,airyfitimg = psf_fitter.psffit(test_airy,returnfitimage=True,vheight=False)
    gaussfitp,gaussfitimg = psf_fitter.psffit(test_airy,psffunction=psf_fitter.twodgaussian,returnfitimage=True,vheight=False)
    print "When fitting, the Gaussian FWHM will be smaller than the Airy FWHM by %f  (inverted, this is %f)" % (gaussfitp[4]*gaussian_fwhm/(airyfitp[4]*airy_fwhm),airyfitp[4]*airy_fwhm/(gaussfitp[4]*gaussian_fwhm))
    print ""

    print "Fitted Airy parameters:     amp %10f xcen %10f ycen %10f width %10f fwhm %10f.  Residual: %10g  Res^2: %10g" % (tuple(airyfitp[1:].tolist())+
            (airyfitp[4]*airy_fwhm,(test_airy-airyfitimg).sum(),((test_airy-airyfitimg)**2).sum()))
    print "Fitted Gaussian parameters: amp %10f xcen %10f ycen %10f width %10f fwhm %10f.  Residual: %10g  Res^2: %10g" % (tuple(gaussfitp[1:].tolist())+
            (gaussfitp[4]*gaussian_fwhm,(test_airy-gaussfitimg).sum(),((test_airy-gaussfitimg)**2).sum()))


    figure(3)
    clf()
    subplot(121)
    imshow((test_airy-airyfitimg), norm=asinh_norm.AsinhNorm())
    colorbar()
    subplot(122)
    imshow((test_airy-gaussfitimg), norm=asinh_norm.AsinhNorm())
    colorbar()

    for ii,noisescale in enumerate(logspace(-2,-0.5,5)):
        gnoise = randn(100,100)*noisescale
        airyfitp,airyfitimg = psf_fitter.psffit(test_airy+gnoise,returnfitimage=True,vheight=False)
        gaussfitp,gaussfitimg = psf_fitter.psffit(test_airy+gnoise,psffunction=psf_fitter.twodgaussian,returnfitimage=True,vheight=False)
        ratio = airyfitp[4]*airy_fwhm/(gaussfitp[4]*gaussian_fwhm)
        print "Fitted Airy parameters:     amp %10f xcen %10f ycen %10f width %10f fwhm %10f ratio %10f.  Residual: %10g  Res^2: %10g" % (tuple(airyfitp[1:].tolist())+(airyfitp[4]*airy_fwhm,ratio,(test_airy-airyfitimg-gnoise).sum(),((test_airy-airyfitimg-gnoise)**2).sum()))
        print "Fitted Gaussian parameters: amp %10f xcen %10f ycen %10f width %10f fwhm %10f ratio %10f.  Residual: %10g  Res^2: %10g" % (tuple(gaussfitp[1:].tolist())+(gaussfitp[4]*gaussian_fwhm,ratio,(test_airy-gaussfitimg-gnoise).sum(),((test_airy-gaussfitimg-gnoise)**2).sum()))
        figure(4+ii)
        clf()
        title("Noise: %0.2g" % noisescale)
        subplot(131)
        title("Noise: %0.2g" % noisescale)
        imshow(test_airy+gnoise,norm=asinh_norm.AsinhNorm())
        colorbar()
        subplot(132)
        imshow(test_airy+gnoise-airyfitimg,norm=asinh_norm.AsinhNorm())
        colorbar()
        subplot(133)
        imshow(test_airy+gnoise-gaussfitimg,norm=asinh_norm.AsinhNorm())
        colorbar()
