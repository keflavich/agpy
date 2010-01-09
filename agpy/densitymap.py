"""
Build a density map out of a list of coordinates
"""

import pyfits
import pyregion
import numpy
from agpy import gaussfitter

def dmregion(header,region,outfits=None,smoothpix=2,clobber=True):
    """
    Given a valid .FITS header with WCS coordinates
    and a ds9 region file, creates a density map of 
    ds9 objects (doesn't filter by shape type).
    Default smoothing is by a sigma=2 pixels gaussian.
    Specify outfits to write to a fits file, otherwise
    returns the map
    """

    x,y = reg_to_xy(region,header)

    return densitymap(header,x,y,outfits=outfits,smoothpix=smoothpix,clobber=clobber)

def reg_to_xy(region,header):
    if region[-3:] == 'reg':
        regfile = pyregion.open(region)
        pix_reg = regfile.as_imagecoord(header)
        x,y = numpy.transpose([rim.coord_list[:2] for rim in pix_reg])
    else:
        print "don't do this it's tough."
        """
        regfile = readcol(region,asdict=True)
        if regfile.has_key('ra'):
            ra = regfile['ra']
            dec = regfile['dec']
            coordsys = 'radec'
        elif regfile.has_key('')
        wcs = pywcs.WCS(header=header)
        """
    return x,y


def densitymap(header,xi,yi,smoothpix=1,outfits=None,clobber=True):
    """
    Generates a source-density map given a region file
    or a list of coordinates
    """


    ny,nx = header['NAXIS2'],header['NAXIS1']
    blankim = numpy.zeros((ny,nx))

    xf = numpy.floor(xi).astype('int')
    xc = numpy.ceil(xi).astype('int')
    yf = numpy.floor(yi).astype('int')
    yc = numpy.ceil(yi).astype('int')
    weight1 = (xi-xf) * (yi-yf)
    weight2 = (xi-xf) * (yc-yi)
    weight3 = (xc-xi) * (yi-yf)
    weight4 = (xc-xi) * (yc-yi)
    blankim[yc,xc] += weight4
    blankim[yc,xf] += weight3
    blankim[yf,xc] += weight2
    blankim[yf,xf] += weight1

    if smoothpix > 1:
        xax,yax = numpy.indices(blankim.shape)
        kernel = gaussfitter.twodgaussian([1,nx/2,ny/2,smoothpix],circle=1,rotate=0,vheight=0)(xax,yax)
        kernelfft = numpy.fft.fft2(kernel)
        imfft = numpy.fft.fft2(blankim)
        dm = numpy.fft.fftshift(numpy.fft.ifft2(kernelfft*imfft).real)
    else:
        dm = blankim

    if outfits:
        hdu = pyfits.PrimaryHDU(dm,header)
        hdu.writeto(outfits,clobber=clobber)

    #for xi,yi in zip(x,y):
    #    
    #    blankim += gaussfitter.twodgaussian([1,xi,yi,smoothpix],circle=1,rotate=0,vheight=0)(xax,yax)

    return dm
