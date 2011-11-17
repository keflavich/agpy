"""
Various codes to work with the initial mass function
"""
import numpy as np

# three codes for dn/dlog(m)
def salpeter(m,alpha=2.35):
    """
    the Salpeter 1955 IMF: dn/dm ~ m^-2.35
    """
    return m**-2.35

def kroupa(m):
    """
    Kroupa 2001 IMF (http://arxiv.org/abs/astro-ph/0009005, http://adsabs.harvard.edu/abs/2001MNRAS.322..231K)
    """
    zeta = (m**0.3 / 0.08**0.3 * 0.08**-1.3)*(m<0.08)
    zeta += (m**-1.3) * (m>=0.08) * (m<0.5)
    zeta += (m**-2.3 / 0.5**-2.3 * 0.5**-1.3) * (m>=0.5)
    return zeta

def chabrier(m):
    """
    Chabrier 2003 IMF
    http://adsabs.harvard.edu/abs/2003PASP..115..763C
    (only valid for m < 1 msun)

    not sure which of these to use...
    """
    # This system MF can be parameterized by the same type of lognormal form as
    # the single MF (eq. [17]), with the same normalization at 1 M⊙, with the
    # coefficients (Chabrier 2003)
    return 0.86 * np.exp(-1*(log(m)-log(0.22))**2/(2*0.57**2))
    # This analytic form for the disk MF for single objects below 1 M⊙, within these uncertainties, is given by the following lognormal form (Chabrier 2003):
    return 0.158 * np.exp(-1*(log(m)-log(0.08))**2/(2*0.69**2))

def schechter(m,A=1,beta=2,m0=100):
    """
    A Schechter function with arbitrary defaults
    """
    return A*m**-beta * np.exp(-m/m0)

#def schechter_inv(m): 
#    """
#    Return p(m)
#    """
#    return scipy.interpolate.interp1d(shfun,arange(.1,20,.01),bounds_error=False,fill_value=20.)

def integrate(fn=kroupa, bins=np.logspace(-2,2,500)):
    xax = (bins[:-1]+bins[1:])/2.
    integral = (bins[1:]-bins[:-1]) * (fn(bins[:-1])+fn(bins[1:])) / 2.

    return xax,integral

def m_integrate(fn=kroupa, bins=np.logspace(-2,2,500)):
    xax = (bins[:-1]+bins[1:])/2.
    integral = xax*(bins[1:]-bins[:-1]) * (fn(bins[:-1])+fn(bins[1:])) / 2.

    return xax,integral

def cumint(fn=kroupa, bins=np.logspace(-2,2,500)):
    xax,integral = integrate(fn,bins)
    return integral.cumsum() / integral.sum()

def m_cumint(fn=kroupa, bins=np.logspace(-2,2,500)):
    xax,integral = m_integrate(fn,bins)
    return integral.cumsum() / integral.sum()
