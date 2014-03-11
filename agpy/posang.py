from astropy import coordinates
from numpy import pi,arctan2,sin,cos,tan

def posang(l1,b1,l2,b2,system='galactic',units='degrees',**kwargs):
    """
    Return the position angle between two points assuming a rectilinear
    coordinate system (I think; at the very least I am making no corrections
    for wcs).

    INPUT:
    longitude1, latitude1, longitude2, latitude2

    Defaults to GALACTIC coordinates.  **kwargs are passed to coords.Position
    """

    if system.lower() == 'galactic':
        pos1 = coordinates.Galactic(l1,b1,unit=('deg','deg'))
        pos2 = coordinates.Galactic(l2,b2,unit=('deg','deg'))
    elif system.lower() in ('radec','fk5','icrs'):
        pos1 = coordinates.ICRS(l1,b1,unit=('deg','deg'))
        pos2 = coordinates.ICRS(l2,b2,unit=('deg','deg'))

    ra1,dec1 = pos1.icrs.ra.deg,pos1.icrs.dec.deg
    ra2,dec2 = pos2.icrs.ra.deg,pos2.icrs.dec.deg

    radiff  = (ra1-ra2)/180.*pi

    angle = arctan2( sin(radiff) , cos(dec1*pi/180.)*tan(dec2*pi/180.) - sin(dec1*pi/180.)*cos(radiff) ) 

    if units == 'degrees':
        return angle/pi*180
    elif units == 'radians':
        return angle
    else:
        raise ValueError("Invalid units: %s" % units)
