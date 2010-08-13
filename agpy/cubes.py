from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array,mean,sum
from math import acos,atan2,tan
from copy import copy
try:
    import pywcs, coords
except ImportError:
    print "cubes.py requires pywcs and coords"

dtor = pi/180.0

def extract_aperture(cube,ap,r_mask=False,wcs=None,coordsys='galactic',wunit='arcsec'):
    """
    Extract an aperture from a data cube.  E.g. to acquire a spectrum
    of an outflow that is extended.

    Cube should have shape [z,y,x], e.g. 
    cube = pyfits.getdata('datacube.fits')

    Apertures are specified in PIXEL units with an origin of 0,0 
    (NOT the 1,1 fits standard!) unless wcs and coordsys are specified
    
    INPUTS:
        wcs - a pywcs.WCS instance associated with the data cube
        coordsys - the coordinate system the aperture is specified in.
            Options are 'celestial' and 'galactic'.  Default is 'galactic'
        wunit - units of width/height.  default 'arcsec', options 'arcmin' and 'degree'

    For a circular aperture, len(ap)=3:
        ap = [xcen,ycen,radius]

    For an elliptical aperture, len(ap)=5:
        ap = [xcen,ycen,height,width,PA]

    Optional inputs:
        r_mask - return mask in addition to spectrum (for error checking?)
    """

    if wcs is not None and coordsys is not None:
        ap = aper_world2pix(ap,wcs,coordsys=coordsys,wunit=wunit)

    if len(ap) == 2:
        sh = cube.shape
        yind,xind = indices(sh[1:3]) # recall that python indices are backwards
        dis = sqrt((xind-ap[0])**2+(yind-ap[1])**2)
        dis3d = repeat(dis[newaxis,:,:],sh[0],axis=0)
        spec = (cube* (dis3d < ap[2])).sum(axis=2).sum(axis=1)
    elif len(ap) == 5:
        yinds,xinds = indices(cube.shape[1:3])
        th = (ap[4])*dtor
        xindr = (xinds-ap[0])*cos(th)  + (yinds-ap[1])*sin(th)
        yindr = (xinds-ap[0])*-sin(th) + (yinds-ap[1])*cos(th)
        ratio = max(ap[2:4])/min(ap[2:4])
        mask = sqrt( (xindr*ratio)**2 + yindr**2) < max(ap[2:4])
        npixinmask = mask.sum()
        mask3d = repeat(mask[newaxis,:,:],cube.shape[0],axis=0)
        spec = (mask3d*cube).sum(axis=2).sum(axis=1) / npixinmask
    else:
        raise Exception("Wrong number of parameters.  Need either 3 parameters "+
                "for a circular aperture or 5 parameters for an elliptical "+ 
                "aperture.")

    if r_mask:
        return spec,mask
    else:
        return spec

def subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,header=None,average=mean):
    """
    Returns a sub-image from a data cube integrated over the specified velocity range

    All units assumed to be pixel units

    cube has dimensions (velocity, y, x)

    """

    xlo = max([xcen-xwidth,0])
    ylo = max([ycen-ywidth,0])
    subim = average(cube[vrange[0]:vrange[1],ylo:ycen+ywidth,xlo:xcen+xwidth],axis=0)

    if header is None:
        return subim

    else:

        hd = header.copy()

        # Header cleanup: must make output 2D.
        try: del hd['CDELT3']
        except: pass
        try: del hd['CD3_3']
        except: pass
        try: del hd['CRVAL3']
        except: pass
        try: del hd['CRPIX3']
        except: pass
        try: del hd['CTYPE3']
        except: pass
        try: del hd['NAXIS3']
        except: pass
        try: del hd['LBOUND3']
        except: pass
        try: del hd['CUNIT3']
        except: pass
        hd['NAXIS']=2

        wcs = pywcs.WCS(header=hd)

        crv1,crv2 = wcs.wcs_pix2sky(xlo,ylo,0)

        hd['CRVAL1'] = crv1[0]
        hd['CRVAL2'] = crv2[0]
        hd['CRPIX1'] = 1
        hd['CRPIX2'] = 1

        return subim,hd

def aper_world2pix(ap,wcs,coordsys='galactic',wunit='arcsec'):
    """
    Converts an elliptical aperture (x,y,width,height,PA) from
    WCS to pixel coordinates given an input wcs (an instance
    of the pywcs.WCS class).  Must be a 2D WCS header.


    """
    convopt = {'arcsec':3600,'arcmin':60,'degree':1}
    try:
        conv = convopt[wunit]
    except:
        raise Exception("Must specify wunit='arcsec','arcmin', or 'degree'")

    if len(wcs.wcs.cdelt) != 2:
        raise Exception("WCS header is not strictly 2-dimensional.  Look for 3D keywords.")
    pos = coords.Position((ap[0],ap[1]),system=coordsys)
    if wcs.wcs.ctype[0][:2] == 'RA':
        ra,dec = pos.j2000()
        corrfactor = cos(dec*dtor)
    elif wcs.wcs.ctype[0][:4] == 'GLON':
        ra,dec = pos.galactic()
        corrfactor=1
    # workaround for a broken wcs.wcs_sky2pix
    try:
        radif = (wcs.wcs.crval[0]-ra)*dtor
        gamma = acos(cos(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)*cos(radif)+sin(dec*dtor)*sin(wcs.wcs.crval[1]*dtor)) / dtor
        theta = atan2( sin(radif) , ( tan(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)-sin(wcs.wcs.crval[1]*dtor)*cos(radif) ) )
        x = -gamma * sin(theta) / wcs.wcs.cd[0,0] + wcs.wcs.crpix[0]
        y = gamma * cos(theta) / wcs.wcs.cd[1,1] + wcs.wcs.crpix[1]
    except:
        radif = (wcs.wcs.crval[0]-ra)*dtor
        gamma = acos(cos(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)*cos(radif)+sin(dec*dtor)*sin(wcs.wcs.crval[1]*dtor)) / dtor
        theta = atan2( sin(radif) , ( tan(dec*dtor)*cos(wcs.wcs.crval[1]*dtor)-sin(wcs.wcs.crval[1]*dtor)*cos(radif) ) )
        x = -gamma * sin(theta) / wcs.wcs.cdelt[0] + wcs.wcs.crpix[0]
        y = gamma * cos(theta) / wcs.wcs.cdelt[1] + wcs.wcs.crpix[1]
    

    #x,y = wcs.wcs_sky2pix(ra,dec,0)  # convert WCS coordinate to pixel coordinate (0 is origin, do not use fits convention)
    try:
        x=x[0]
        y=y[0]
    except:
        pass
    # cd is default, cdelt is backup
    try:
        width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
        height = ap[3] / conv / abs(wcs.wcs.cd[0,0])
    except:
        width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
        height = ap[3] / conv / abs(wcs.wcs.cdelt[0])
    PA = ap[4] 
    apold = copy(ap)
    ap = [x,y,width,height,PA]

    return ap





