from numpy import sqrt,repeat,indices,newaxis,pi,cos,sin,array,mean,sum,nansum
from math import acos,atan2,tan
import copy
import pyfits
import pyregion
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

def integ(file,vrange,xcen=None,xwidth=None,ycen=None,ywidth=None,**kwargs):
    """
    wrapper of subimage_integ that defaults to using the full image
    """
    if isinstance(file,pyfits.PrimaryHDU):
        header = file.header
        cube = file.data
    elif isinstance(file,pyfits.HDUList):
        header = file[0].header
        cube = file[0].data
    else:
        file = pyfits.open(file)
        header = file[0].header
        cube = file[0].data

    if None in [xcen,xwidth,ycen,ywidth]:
        xcen = header['NAXIS1'] / 2
        xwidth = xcen
        ycen = header['NAXIS2'] / 2
        ywidth = ycen

    return subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,**kwargs)

def subimage_integ(cube,xcen,xwidth,ycen,ywidth,vrange,header=None,average=mean,dvmult=False,units="pixels"):
    """
    Returns a sub-image from a data cube integrated over the specified velocity range

    All units assumed to be pixel units

    cube has dimensions (velocity, y, x)

    xwidth and ywidth are "radius" values, i.e. half the length that will be extracted

    if dvmult is set, multiple the average by DV (this is useful if you set
    average=sum and dvmul=True to get an integrated value)

    """

    if header:
        flathead = flatten_header(header.copy())
        wcs = pywcs.WCS(header=flathead)
        if header.get('CD3_3'): CD3 = header.get('CD3_3')
        else: CD3 = header.get('CDELT3')

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
        zrange = ( array(vrange)-header.get('CRVAL3') ) / CD3 - 1 + header.get('CRPIX3')
    else:
        print "Can only use wcs if you pass a header."

    subim = average(cube[zrange[0]:zrange[1],ylo:yhi,xlo:xhi],axis=0)
    if dvmult and CD3: subim *= CD3
    elif dvmult: print "Error: could not multiply by dv; CD3=",CD3

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


def getspec(lon,lat,rad,cube,header,r_fits=True,inherit=True,wunit='arcsec'):
    """
    Given a longitude, latitude, aperture radius (arcsec), and a cube file,
    return a .fits file or a spectrum.
    
    lon,lat - longitude and latitude center of a circular aperture in WCS coordinates
    rad     - radius (default degrees) of aperture
    """

    convopt = {'arcsec':1.0,'arcmin':60.0,'degree':3600.0}

    flathead = flatten_header(header)
    wcs = pywcs.WCS(flathead)
    if wcs.wcs.ctype[0][:2] == 'RA':
      coordsys='celestial'
    elif wcs.wcs.ctype[0][:4] == 'GLON':
      coordsys='galactic'
    spec = extract_aperture(cube,[lon,lat,rad],wcs=wcs,
            coordsys=coordsys,wunit=wunit)

    if nansum(spec) == 0:
        print "Total of extracted spectrum was zero. lon,lat,rad: ",lon,lat,rad #  Tracing to find your problem."
        #import pdb; pdb.set_trace()

    if r_fits:
        if inherit:
            newhead = header.copy()
        else:
            newhead = pyfits.Header()
        try:
            newhead.update('CD1_1',header['CD3_3'])
        except KeyError:
            newhead.update('CD1_1',header['CDELT3'])
        newhead.update('CRPIX1',header['CRPIX3'])
        newhead.update('CRVAL1',header['CRVAL3'])
        try:
            newhead.update('CTYPE1',header['CTYPE3'])
        except KeyError:
            newhead.update('CTYPE1',"VRAD")
        try:
            newhead.update('CUNIT1',header['CUNIT3'])
        except KeyError:
            print "Header did not contain CUNIT3 keyword.  Defaulting to km/s"
            newhead.update('CUNIT1',"km/s")
        newhead.update('BUNIT',header['BUNIT'])
        newhead.update('APGLON',lon)
        newhead.update('APGLAT',lat)
        newhead.update('APRAD',rad*convopt[wunit],comment='arcseconds') # radius in arcsec
        newfile = pyfits.PrimaryHDU(data=spec,header=newhead)
        return newfile
    else:
        return spec

def getspec_reg(cubefilename,region):
    """
    Aperture extraction from a cube using a pyregion circle region

    The region must be in the same coordinate system as the cube header
    """

    ds9tocoords = {'fk5':'celestial','galactic':'galactic','icrs':'celestial'}

    if region.name != 'circle':
        raise Exception("Only circular apertures are implemented so far")

    l,b,r = region.coord_list
    #pos = coords.Position([l,b],system=ds9tocoords[region.coord_format])
    cubefile = pyfits.open(cubefilename)
    header = cubefile[0].header
    cube = cubefile[0].data
    if len(cube.shape) == 4: cube = cube[0,:,:,:]

    sp = getspec(l,b,r,cube,header,wunit='degree')

    return sp

def smooth_cube(cube,cubedim=0,parallel=True,**kwargs):
    """
    parallel-map the smooth function
    """
    from convolve import smooth
    from contributed import parallel_map

    if cubedim != 0:
        cube = cube.swapaxes(0,cubedim)

    cubelist = [cube[ii,:,:] for ii in xrange(cube.shape[0])]

    Psmooth = lambda C: smooth(C,**kwargs)

    if parallel:
        smoothcube = array(parallel_map(Psmooth,cubelist))
    else:
        smoothcube = array(map(Psmooth,cubelist))
    
    if cubedim != 0:
        smoothcube = smoothcube.swapaxes(0,cubedim)

    return smoothcube


