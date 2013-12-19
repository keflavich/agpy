import numpy as np
from numpy import pi

pc = 3.086e16 # m  (RATRAN uses meters)
au = 1.496e13 # cm
msun = 2e33 # grams
mh = 1.67e-24
G = 6.67e-8 # CGS

def mt03(r, column=1e22, krho=1.5, mass=1.0, rhoouter=1e3, zh2=2.8,mh=1.66053886e-24, r_inner=au):
    """
    Remember, radius is given in METERS
    r_inner is given in cm
    """
    r = r*1e2 # convert to cm
    r[r<r_inner] = r_inner
    r_outer = (mass*msun/(pi * column * mh * zh2))**0.5
    A = (4-krho)*mass*msun / (4/3. * pi * r_outer**(4-krho)) 
    # equivalent to prev. line A = mass*msun / (r_outer**(4-krho)) * (3*(4-krho)/(4.*pi))

    rho = ((A * r**-krho / (zh2*mh) ) * (r<r_outer) +
            rhoouter*zh2*mh * (r>r_outer))

    nh2 = rho / zh2 / mh

    return nh2


def plummer(r,mtot=0,a=0,zh2=2.8,mh=1.66053886e-24,ismeters=True):
    """
    Return the density given a Plummer profile
    """
    
    scale = 1e-6 if ismeters else 1.0
    rho = 3 * mtot / (4*np.pi*a**3) * ( 1. + r**2/a**2 )**(-2.5)
    nh2 = rho / (mh*zh2)
    return nh2 * scale

def broken_powerlaw(r,rbreak=1.0,power=-2.0,n0=1e5):
    """
    Return the density of a broken power law density profile where
    the central density is flat
    """

    if np.isscalar(r):
        if r < rbreak:
            return n0
        else:
            # n(rbreak) = n0
            n = n0 * (r/rbreak)**power
            return n
    else:
        narr = np.zeros(r.shape)
        narr[r<rbreak] = n0
        # this was previously
        # narr[r>=rbreak] = n0 * (r/rbreak)**power
        # which is wrong on two levels - first, the math is wrong, but second,
        # the assignment is supposed to be disallowed because of shape mismatch
        # numpy has failed me.
        narr[r>=rbreak] = n0 * (r[r>=rbreak]/rbreak)**power
        return narr

def bonnorebert(r,ncenter=1e5,ximax=6.9,viso=0.24,zh2=2.8):
    """
    approximation to a bonnor-ebert sphere using broken power law
    6.9 taken from Alves, Lada, Lada B68 Nature paper 2001
    viso is for zh2=2.8, T=20K
    """
    # cm
    rbreak = ximax*(viso*1e5) / np.sqrt(4*np.pi*G*ncenter*mh*zh2)
    # to m 
    rbreak /= 100.
    print "bonnorebert rbreak: %g " % rbreak

    return broken_powerlaw(r,rbreak=rbreak,n0=ncenter)

def king(r, rc, rt, sigma_0, alpha=2):
    """
    See http://iopscience.iop.org/1538-3881/139/6/2097/fulltext/

    Parameters
    ----------
    r: float
        radius
    rc: float
        core radius
    rt: float
        truncation radius
    sigma_0: float
        central density
    """

    def z(x):
        return 1/(1+(x/rc)**2)**(1./alpha)

    term1 = (1 - z(rt))**-alpha
    term2 = (z(r) - z(rt))**alpha
    sigma = sigma_0 * term1 * term2

    return sigma
