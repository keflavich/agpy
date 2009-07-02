from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array
from math import acos,atan2,tan
from copy import copy
import pywcs, coords

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
        [NOTE: DOES NOT USE PYWCS TO CALCULATE PIXEL COORDINATES!  They crash.
        XY computation pilfered from gcirc.pro and posang.pro in the astron idl library]

    For a circular aperture, len(ap)=3:
        ap = [xcen,ycen,radius]

    For an elliptical aperture, len(ap)=5:
        ap = [xcen,ycen,height,width,PA]

    Optional inputs:
        r_mask - return mask in addition to spectrum (for error checking?)
    """
    convopt = {'arcsec':3600,'arcmin':60,'degree':1}
    try:
        conv = convopt[wunit]
    except:
        raise Exception("Must specify wunit='arcsec','arcmin', or 'degree'")

    if wcs is not None and coordsys is not None:
        pos = coords.Position((ap[0],ap[1]),system=coordsys)
        if wcs.wcs.ctype[0][:2] == 'RA':
            ra,dec = pos.j2000()
            corrfactor = cos(dec*dtor)
        elif wcs.wcs.ctype[0][:4] == 'GLON':
            ra,dec = pos.galactic()
            corrfactor=1
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

        #x,y = wcs.wcs_sky2pix(array([ra]),array([dec]),0)  # convert WCS coordinate to pixel coordinate (0 is origin, do not use fits convention)
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cd[0,0])
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cdelt[0])
        PA = ap[4] 
        apold = copy(ap)
        ap = [x,y,width,height,PA]

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
        print "Extracted ellipse"
    else:
        raise Exception("Wrong number of parameters.  Need either 3 parameters "+
                "for a circular aperture or 5 parameters for an elliptical "+ 
                "aperture.")

    if r_mask:
        return spec,mask
    else:
        return spec

    
    


