from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array,mean,sum,nansum
from math import acos,atan2,tan
import copy
import pyfits
try:
    import pywcs, coords
except ImportError:
    print "cubes.py requires pywcs and coords"

dtor = pi/180.0

def flatten_header(header):
    """
    Attempt to turn an N-dimensional fits header into a 2-dimensional header
    Turns all CRPIX[>2] etc. into new keywords with suffix 'A'

    header must be a pyfits.Header instance
    """

    if not isinstance(header,pyfits.Header):
        raise Exception("flatten_header requires a pyfits.Header instance")

    newheader = header.copy()

    for key in newheader.keys():
        try:
            if int(key[-1]) >= 3 and key[:2] in ['CD','CR','CT','CU','NA']:
                newheader.rename_key(key,'A'+key,force=True)
        except ValueError:
            # if key[-1] is not an int
            pass
        except IndexError:
            # if len(key) < 2
            pass
    newheader.update('NAXIS',2)

    return newheader

def extract_aperture(cube,ap,r_mask=False,wcs=None,coordsys='galactic',wunit='arcsec'):
    """
    Extract an aperture from a data cube.  E.g. to acquire a spectrum
    of an outflow that is extended.

    Cube should have shape [z,y,x], e.g. 
    cube = pyfits.getdata('datacube.fits')

    Apertures are specified in PIXEL units with an origin of 0,0 (NOT the 1,1
    fits standard!) unless wcs and coordsys are specified
    
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

    if len(ap) == 3:
        sh = cube.shape
        yind,xind = indices(sh[1:3]) # recall that python indices are backwards
        dis = sqrt((xind-ap[0])**2+(yind-ap[1])**2)
        mask = dis < ap[2]
    elif len(ap) == 5:
        yinds,xinds = indices(cube.shape[1:3])
        th = (ap[4])*dtor
        xindr = (xinds-ap[0])*cos(th)  + (yinds-ap[1])*sin(th)
        yindr = (xinds-ap[0])*-sin(th) + (yinds-ap[1])*cos(th)
        ratio = max(ap[2:4])/min(ap[2:4])
        mask = sqrt( (xindr*ratio)**2 + yindr**2) < max(ap[2:4])
    else:
        raise Exception("Wrong number of parameters.  Need either 3 parameters "+
                "for a circular aperture or 5 parameters for an elliptical "+ 
                "aperture.")

    npixinmask = mask.sum()
    mask3d = repeat(mask[newaxis,:,:],cube.shape[0],axis=0)
    spec = nansum(nansum((cube*mask3d),axis=2),axis=1) / npixinmask

    if r_mask:
        return spec,mask
    else:
        return spec

def subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,header=None,average=mean,units="pixels"):
    """
    Returns a sub-image from a data cube integrated over the specified velocity range

    All units assumed to be pixel units

    cube has dimensions (velocity, y, x)

    xwidth and ywidth are "radius" values, i.e. half the length that will be extracted

    """

    if header:
        flathead = flatten_header(header.copy())
        wcs = pywcs.WCS(header=flathead)

    if units=="pixels":
        xlo = int( max([xcen-xwidth,0])              )
        ylo = int( max([ycen-ywidth,0])              )
        xhi = int( min([xcen+xwidth,cube.shape[2]])  )
        yhi = int( min([ycen+ywidth,cube.shape[1]])  )
        zrange = vrange
    elif units=="wcs" and header:
        newxcen,newycen = wcs.wcs_sky2pix(xcen,ycen,0)
        try:
            newxwid,newywid = xwidth / abs(wcs.wcs.cd[0,0]), ywidth / abs(wcs.wcs.cd[1,1])
        except AttributeError:
            newxwid,newywid = xwidth / abs(wcs.wcs.cdelt[0]), ywidth / abs(wcs.wcs.cdelt[1])
        xlo = int( max([newxcen-newxwid,0]) )
        ylo = int( max([newycen-newywid,0]) )
        xhi = int( min([newxcen+newxwid,cube.shape[2]]) )
        yhi = int( min([newycen+newywid,cube.shape[1]]) )
        if header.get('CD3_3'):
            zrange = ( array(vrange)-header.get('CRVAL3') ) / header.get('CD3_3') - 1 + header.get('CRPIX3')
        else:
            zrange = ( array(vrange)-header.get('CRVAL3') ) / header.get('CDELT3') - 1 + header.get('CRPIX3')
    else:
        print "Can only use wcs if you pass a header."

    subim = average(cube[zrange[0]:zrange[1],ylo:yhi,xlo:xhi],axis=0)

    if header is None:
        return subim
    else:
        crv1,crv2 = wcs.wcs_pix2sky(xlo,ylo,0)

        flathead['CRVAL1'] = crv1[0]
        flathead['CRVAL2'] = crv2[0]
        flathead['CRPIX1'] = 1
        flathead['CRPIX2'] = 1

        return subim,flathead

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
    if len(ap) == 5:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cd[0,0])
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
            height = ap[3] / conv / abs(wcs.wcs.cdelt[0])
        PA = ap[4] 
        apold = copy.copy(ap)
        ap = [x,y,width,height,PA]
    elif len(ap) == 3:
        try:
            width  = ap[2] / conv / abs(wcs.wcs.cd[0,0])  # first is width, second is height in DS9 PA convention
        except:
            width  = ap[2] / conv / abs(wcs.wcs.cdelt[0])  # first is width, second is height in DS9 PA convention
        apold = copy.copy(ap)
        ap = [x,y,width]

    return ap





