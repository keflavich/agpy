from numpy import sign

def ratos(ra):
    r = ra/15.
    deg = r - (r % 1)
    min = (r % 1) * 60.
    sec = (min % 1) * 60
    rs = "%02.0f:%02.0f:%05.2f" % (deg,min,sec)
    return rs

def dectos(dec):
    deg = sign(dec) * (abs(dec) - (abs(dec) % 1))
    min = (abs(dec) % 1) * 60.
    sec = (min % 1) * 60
    decs = "%+02.0f:%02.0f:%05.2f" % (deg,min,sec)
    return decs

