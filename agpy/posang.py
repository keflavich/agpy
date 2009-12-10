import coords
from numpy import pi,arctan2,sin,cos,tan

def posang(l1,b1,l2,b2,system='galactic'):
    pos1 = coords.Position([l1,b1],system=system)
    ra1,dec1 = pos1.j2000()
    pos2 = coords.Position([l2,b2],system=system)
    ra2,dec2 = pos2.j2000()

    radiff  = (ra1-ra2)/180.*pi

    angle = arctan2( sin(radiff) , cos(dec1*pi/180.)*tan(dec2*pi/180.) - sin(dec1*pi/180.)*cos(radiff) ) /pi*180

    return angle
