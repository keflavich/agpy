"""
Integrate a cube over some velocities and plot it with aplpy.  Most of this
functionality is subsumed by :mod:cubes
"""
try:
    import aplpy
except ImportError:
    print "aplpy is required for collapseplot"
import pyfits
from pylab import *
for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'

def collapseplot(filename, vrange=[0,10], outfilename=None, contourfile=None,
        xrange=None, yrange=None, imfile=None, clobber=True, swapcontour=False, 
        hdu=0, vmin=None, vmax=None, integrate=False, **kwargs):
    """
    Collapses a data cube over velocities specified by vrange
    Inputs:
    filename - name of data cube (.fits)
    vrange - velocity range (units depend on header)

    Optional inputs:
    outfilename - will write collapsed cube to an outfile.
        clobber - if outfilename set, allows clobbering
    contourfile - overlays contours from this file.  Accepts **kwargs
    swapcontour - switch which image will be greyscale and which will be contours
    imfile - saves file as an image
    xrange, yrange - crops image to [xmin,xmax] and [ymin,ymax] where x,y
        are in the header coordinates
    vmin,vmax - grayscale for plotting
    """

    file = collapsecube(filename,vrange,outfilename=outfilename,clobber=clobber,integrate=integrate,hdu=hdu)

    if xrange or yrange:
        file = cropfits(file,xrange=xrange,yrange=yrange)

    
    if contourfile:
        cf = pyfits.open(contourfile)
        cf = cropfits(cf,xrange=xrange,yrange=yrange)
        if swapcontour:
            fig = aplpy.FITSFigure(cf)
            fig.show_grayscale(vmin=vmin,vmax=vmax)
            fig.show_contour(file, **kwargs)
    else:
        fig = aplpy.FITSFigure(file)
        fig.show_grayscale(vmin=vmin,vmax=vmax)
        if contourfile: fig.show_contour(cf, **kwargs)

    if imfile:
        fig.save(imfile,dpi=150)

    show()

    return fig
    
    
def cropfits(fitsfile,xrange=None,yrange=None,hdu=0):
    """
    Crop a fits file in the range specified by xrange,yrange 
    in header units
    """
    
    im = fitsfile[hdu].data
    header = fitsfile[hdu].header 

    try:
        dx,x0,x0pix = header['CD1_1'],header['CRVAL1'],header['CRPIX1']
        dy,y0,y0pix = header['CD2_2'],header['CRVAL2'],header['CRPIX2']
    except:
        try: 
            dx,x0,x0pix = header['CDELT1'],header['CRVAL1'],header['CRPIX1']
            dy,y0,y0pix = header['CDELT2'],header['CRVAL2'],header['CRPIX2']
        except: 
            raise Exception('Failed to read header')

    x = (arange(header['NAXIS1'])-(x0pix-1.))*dx + x0
    y = (arange(header['NAXIS2'])-(y0pix-1.))*dy + y0

    if xrange:
        xmin = argmin(abs(x-xrange[1]))
        xmax = argmin(abs(x-xrange[0]))
        header['CRVAL1'] = x[xmin]
        header['CRPIX1'] = 1.0
        im = im[:,xmin:xmax]
    if yrange:
        ymin = argmin(abs(y-yrange[0]))
        ymax = argmin(abs(y-yrange[1]))
        header['CRVAL2'] = y[ymin]
        header['CRPIX2'] = 1.0
        im = im[ymin:ymax,:]

    fitsfile[0].data = im
    header['NAXIS1'] = im.shape[1]
    header['NAXIS2'] = im.shape[0]

    return fitsfile

def collapsecube(filename,vrange,outfilename=None,clobber=True,integrate=True,hdu=0):
    """ 
    Collapses a datacube over some velocity range and returns the
    pyfits HDU that includes it.  Optionally writes to disk.
    filename - input DATA CUBE file name
    vrange - velocity or frequency range over which to sum

    Optional inputs:
    outfilename - file name to write to (default None)
    clobber - overwrite extant filename (default True)
    integrate - is output an integral (i.e. sum * delta-V) [default True]

    """
    file = pyfits.open(filename)
    cube = file[hdu].data
    header = file[hdu].header

    try:
        dv,v0,v0pix = header['CD3_3'],header['CRVAL3'],header['CRPIX3']
    except:
        try: 
            dv,v0,v0pix = header['CDELT3'],header['CRVAL3'],header['CRPIX3']
            header.update('CD1_1',header['CDELT1'])
            header.update('CD2_2',header['CDELT2'])
        except:
            raise Exception("Failed to read 3rd dimension of data cube.")
    vel = (arange(cube.shape[0])-(v0pix-1.))*dv + v0

    im = cube[((vel>vrange[0])*(vel<vrange[1])),:,:].sum(axis=0)
    if integrate:
        im *= dv

    file[0].data = im

    header['NAXIS'] = 2
    del header['CRVAL3']
    del header['CRPIX3']
    del header['NAXIS3']
    try: del header['CDELT3']
    except: pass
    try: del header['CD3_3']
    except: pass
    
    # Deal with weird header types
    if header['CTYPE1'][-3:] == 'GLS':
        header['CTYPE1'] = header['CTYPE1'][:-3]+'CAR'
    if header['CTYPE2'][-3:] == 'GLS':
        header['CTYPE2'] = header['CTYPE2'][:-3]+'CAR'

    if outfilename:
        file.writeto(outfilename, clobber=clobber)

    return file

