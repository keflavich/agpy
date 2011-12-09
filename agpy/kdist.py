try: 
    import coords
except ImportError:
    print "Can't use kdist without the 'coords' package"
from numpy import sqrt, abs, pi, cos, sin, max, ones, array

def kdist(l, b, vin, near=True,r0=8.4e3,v0=2.54e2,dynamical=False,
        kinematic=True,regular=False,rrgal=False,verbose=False,
        inverse=False,silent=False, returnvtan=False):
    """
     NAME:
       KINDIST 
     PURPOSE:
       To return the distance to an object given l,b,v
    
     CALLING SEQUENCE:
       dist = KDIST (L, B, V)
    
     INPUTS:
       L, B -- Galactic Longitude and Latitude (decimal degrees)
       V - Velocity w.r.t. LSR in km/s
     KEYWORD PARAMETERS:
       /NEAR, /FAR -- Report the near/far kinematic distances for Q1 and
                      Q4 data.
       RO, VO -- Force values for galactocentric distance for sun and
                 velocity of the LSR around the GC.  Default to 8.4 kpc
                 and 254 km/s (Reid et al., 2009)
       RGAL -- Named keyword containing galactocentric radius of sources.
       rrgal  - return galactocentric distance in addition to distance from us
       /DYNAMICAL -- Use the dynamical definition of the LSR
       /KINEMATIC -- Use the kinematic definition of the LSR (default)
       /REGULAR -- Do not apply the rotation correction for High mass
                   star forming regions.
        INVERSE -- If set, pass DISTANCE instead of velocity, and output is
                   velocity
        returnvtan - if set, return the tanent velocity and ignore the input
                    velocity
     OUTPUTS:
       DIST -- the kinematic distance in units of R0 (defaults to pc).
    
     MODIFICATION HISTORY:
    
           Fri Feb 27 00:47:18 2009, Erik <eros@orthanc.local>
    		 Adapted from kindist.pro
                 Translated from IDL to Python by Adam Ginsburg (adam.ginsburg@colorado.edu)

    """

    dtor = pi/180.

    if regular: vs = 0.0 
    else: vs=15.0

    if kinematic or not(dynamical):
        solarmotion_ra = ((18+03/6e1+50.29/3.6e3)*15)
        solarmotion_dec = (30+0/6e1+16.8/3.6e3)
        solarmotion_mag = 20.0
    else:
        solarmotion_ra = ((17+49/6e1+58.667/3.6e3)*15)
        solarmotion_dec = (28+7/6e1+3.96/3.6e3)
        solarmotion_mag = 16.55294

    cg = coords.Position((l,b),system='galactic')
    solarmotion = coords.Position((solarmotion_ra,solarmotion_dec))
    #  ra,dec = cg.j2000()
    #  gcirc, 2, solarmotion_ra, solarmotion_dec, ra, dec, theta
    theta = cg.angsep(solarmotion).arcsec()

    vhelio = vin-solarmotion_mag*cos(theta/206265.)

    # UVW from Dehnen and Binney
    bigu = 10.0
    bigv = 5.23
    bigw = 7.17

    v = vhelio+(bigu*cos(l*dtor)+bigv*sin(l*dtor))*cos(b*dtor)+bigw*sin(b*dtor)

    # Compute tangent distance and velocity
    rtan = r0*(cos(l*dtor))/(cos(b*dtor))
    vTEMP = (1/sin(l*dtor) - v0/(v0-vs)) * ((v0-vs)*sin(l*dtor)*cos(b*dtor))
    vhelioTEMP = vTEMP - ((bigu*cos(l*dtor)+bigv*sin(l*dtor))*cos(b*dtor)+bigw*sin(b*dtor))
    vtan = vhelioTEMP+solarmotion_mag*cos(theta/206265.)
    if returnvtan:
        return vtan

    # This is r/r0
    null = (v0/(v0-vs)+v/((v0-vs)*sin(l*dtor)*cos(b*dtor)))**(-1)

    if inverse:
        radical = cos(l*dtor) - cos(b*dtor) * vin / r0 
        null = sqrt(1 - cos(l*dtor)**2 + radical**2)
        v = (1/null - v0/(v0-vs)) * ((v0-vs)*sin(l*dtor)*cos(b*dtor))
        vhelio = v - ((bigu*cos(l*dtor)+bigv*sin(l*dtor))*cos(b*dtor)+bigw*sin(b*dtor))
        vlsr = vhelio+solarmotion_mag*cos(theta/206265.)
        return vlsr
    else:
        if vin > vtan:
            if not silent:
                print "Velocity is greater than tangent velocity v=%f.  Returning tangent distance." % vtan
            if rrgal: return rtan,null*r0
            return rtan
        #  The > 0 traps things near the tangent point and sets them to the
        #  tangent distance.  So quietly.  Perhaps this should pitch a flag?
        radical = max(sqrt(((cos(l*dtor))**2-(1-null**2)) ),0)

        fardist = r0*(cos(l*dtor)+radical)/(cos(b*dtor))

        neardist = r0*(cos(l*dtor)-radical)/(cos(b*dtor))

    rgal = null*r0
    ind = (abs(l-180) < 90)
    if ind.sum() > 1: neardist[ind] = fardist[ind]
    elif ind==True: neardist = fardist

    if not(near): dist = fardist
    else: dist = neardist

    if verbose:
        print "radical: %f  null: %f  vin: %f  v: %f  vhelio: %f rgal: %f  neardist: %f  fardist: %f" % (radical,null,vin,v,vhelio,rgal,neardist,fardist)

    if rrgal: return abs(dist),abs(rgal)
    return abs(dist)

def vector_kdist(x,y,z,**kwargs):

    if type(z)==type(1) or type(z)==type(1.0):
        z = z*ones(len(x))
    v = []
    for i,j,k in array([x,y,z]).T:
        v.append( kdist(i,j,k,**kwargs) )
    return array(v)

def threekpcarm(longitude,radius=3.0,center_distance=8.5):
    return sqrt(radius**2+center_distance**2-2*radius*center_distance*cos( (90-3*longitude) / 180. * pi ))
