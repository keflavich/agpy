"""
Generate a cutout image from a .fits file
"""
import pyfits
import numpy
import pywcs

def cutout(file,xc,yc,xw=25,yw=25,units='pixels',outfile=None,clobber=True):
    """
    Inputs:
        file  - .fits filename or pyfits HDUList (must be 2D)
        xc,yc - x and y coordinates in the fits files' coordinate system
        xw,yw - x and y width 
        units - specify units to use: either pixels or wcs
        outfile - optional output file
    """

    if isinstance(file,str):
        file = pyfits.open(file)
        opened=True
    elif isinstance(file,pyfits.HDUList):
        opened=False
    else:
        raise Exception("cutout: Input file is wrong type (string or HDUList are acceptable).")

    head = file[0].header.copy()

    if head['NAXIS'] > 2:
        raise Exception("Too many (%i) dimensions!" % head['NAXIS'])
    try:
        cd1 = head['CDELT1']
        cd2 = head['CDELT2']
    except KeyError:
        try:
            cd1 = head['CD1_1']
            cd2 = head['CD2_2']
        except KeyError:
            raise Exception("No CD or CDELT keywords in header")

    lonarr = ((numpy.arange(head['NAXIS1'])-head['CRPIX1'])*cd1 + head['CRVAL1'] )
    latarr = ((numpy.arange(head['NAXIS2'])-head['CRPIX2'])*cd2 + head['CRVAL2'] )

    wcs = pywcs.WCS(head)

    #xx = numpy.argmin(numpy.abs(xc-lonarr))
    #yy = numpy.argmin(numpy.abs(yc-latarr))
    xx,yy = wcs.wcs_sky2pix(xc,yc,0)


    if units=='pixels':
        xmin,xmax = numpy.max([0,xx-xw]),numpy.min([head['NAXIS1'],xx+xw])
        ymin,ymax = numpy.max([0,yy-yw]),numpy.min([head['NAXIS2'],yy+yw])
    elif units=='wcs':
        xmin,xmax = numpy.max([0,xx-xw/numpy.abs(cd1)]),numpy.min([head['NAXIS1'],xx+xw/numpy.abs(cd1)])
        ymin,ymax = numpy.max([0,yy-yw/numpy.abs(cd2)]),numpy.min([head['NAXIS2'],yy+yw/numpy.abs(cd2)])
    else:
        raise Exception("Can't use units %s." % units)

    if xmax < 0 or ymax < 0:
        raise ValueError("Coordinate is outside of map: %f,%f." % (xmax,ymax))

    img = file[0].data[ymin:ymax,xmin:xmax]

    head['CRPIX1']-=xmin
    head['CRPIX2']-=ymin
    head['NAXIS1']=img.shape[1]
    head['NAXIS2']=img.shape[0]

    newfile = pyfits.PrimaryHDU(data=img,header=head)

    if isinstance(outfile,str):
        newfile.writeto(outfile,clobber=clobber)

    if opened:
        file.close()

    return newfile
