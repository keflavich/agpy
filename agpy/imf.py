"""
Various codes to work with the initial mass function
"""
import numpy as np
import types # I use typechecking.  Is there a better way to do this?  (see inverse_imf below)

# three codes for dn/dlog(m)
def salpeter(m,alpha=2.35, integral=False):
    """
    the Salpeter 1955 IMF: dn/dm ~ m^-2.35
    """
    if integral: alpha -= 1
    return m**(-alpha)

def kroupa(m,integral=False):
    """
    Kroupa 2001 IMF (http://arxiv.org/abs/astro-ph/0009005, http://adsabs.harvard.edu/abs/2001MNRAS.322..231K)
    """
    exp1 = 0.3
    exp2 = 1.3
    exp3 = 2.3
    if integral: 
        exp1 += 1
        exp2 -= 1
        exp3 -= 1
    zeta = (m**exp1 / 0.08**exp1 * 0.08**-exp2)*(m<0.08)
    zeta += (m**-exp2) * (m>=0.08) * (m<0.5)
    zeta += (m**-exp3 / 0.5**-exp3 * 0.5**-exp2) * (m>=0.5)
    return zeta

def chabrier(m, integral=False):
    """
    Chabrier 2003 IMF
    http://adsabs.harvard.edu/abs/2003PASP..115..763C
    (only valid for m < 1 msun)

    not sure which of these to use...

    integral is NOT IMPLEMENTED
    """
    if integral: print "Chabrier integral NOT IMPLEMENTED"
    # This system MF can be parameterized by the same type of lognormal form as
    # the single MF (eq. [17]), with the same normalization at 1 Msun, with the
    # coefficients (Chabrier 2003)
    return 0.86 * np.exp(-1*(np.log10(m)-np.log10(0.22))**2/(2*0.57**2))
    # This analytic form for the disk MF for single objects below 1 Msun, within these uncertainties, is given by the following lognormal form (Chabrier 2003):
    return 0.158 * np.exp(-1*(np.log10(m)-np.log10(0.08))**2/(2*0.69**2))

def schechter(m,A=1,beta=2,m0=100, integral=False):
    """
    A Schechter function with arbitrary defaults
    (integral may not be correct - exponent hasn't been dealt with at all)
    """
    if integral: beta -= 1
    return A*m**-beta * np.exp(-m/m0)

try: 
    import scipy
    def schechter_cdf(m,A=1,beta=2,m0=100,mmin=10,mmax=None,npts=1e4):
        """
        Return the CDF value of a given mass for a set mmin,mmax
        mmax will default to 10 m0 if not specified
        
        Analytic integral of the Schechter function:
        http://www.wolframalpha.com/input/?i=integral%28x^-a+exp%28-x%2Fm%29+dx%29
        """
        if mmax is None:
            mmax = 10*m0
        
        # integrate the CDF from the minimum to maximum 
        # undefined posint = -m0 * mmax**-beta * (mmax/m0)**beta * scipy.special.gammainc(1-beta, mmax/m0)
        # undefined negint = -m0 * mmin**-beta * (mmin/m0)**beta * scipy.special.gammainc(1-beta, mmin/m0)
        posint = -mmax**(1-beta) * scipy.special.expn(beta, mmax/m0)
        negint = -mmin**(1-beta) * scipy.special.expn(beta, mmin/m0)
        tot = posint-negint

        # normalize by the integral
        # undefined ret = (-m0 * m**-beta * (m/m0)**beta * scipy.special.gammainc(1-beta, m/m0)) / tot
        ret = (-m**(1-beta) * scipy.special.expn(beta, m/m0) - negint)/ tot

        return ret

    def sh_cdf_func(**kwargs):
        return lambda x: schechter_cdf(x, **kwargs)
except ImportError:
    pass




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

massfunctions = {'kroupa':kroupa, 'salpeter':salpeter, 'chabrier':chabrier, 'schechter':schechter}
# salpeter and schechter selections are arbitrary
mostcommonmass = {'kroupa':0.08, 'salpeter':0.01, 'chabrier':0.23, 'schecter':0.01}

def inverse_imf(p, nbins=1000, mmin=0.03, mmax=120, massfunc='kroupa'):
    """
    Inverse mass function

    massfunc can be 'kroupa', 'chabrier', 'salpeter', 'schechter', or a function
    """
 
    masses = np.logspace(np.log10(mmin),np.log10(mmax),nbins)
    if type(massfunc) is types.FunctionType:
        mf = massfunc(masses, integral=True)
    elif type(massfunc) is str:
        mf = massfunctions[massfunc](masses, integral=True)
    else:
        raise ValueError("massfunc must either be a string in the set %s or a function" % (",".join(massfunctions.keys())))
    mfcum = mf.cumsum()
    mfcum /= mfcum.max() # normalize to sum (cdf)

    return np.interp(p, mfcum, masses)

def make_cluster(mcluster, massfunc='kroupa', verbose=False, silent=False, tolerance=0.5, **kwargs):
    """
    Sample from an IMF to make a cluster.  Returns the masses of all stars in the cluster
    
    massfunc must be a string 
    tolerance is how close the cluster mass must be to the requested mass.  
    If the last star is greater than this tolerance, the total mass will not be within
    tolerance of the requested

    kwargs are passed to `inverse_imf`
    """

    # use most common mass to guess needed number of samples
    nsamp = mcluster / mostcommonmass[massfunc]
    masses = inverse_imf(np.random.random(nsamp), massfunc=massfunc, **kwargs)

    mtot = masses.sum()
    if verbose: print "%i samples yielded a cluster mass of %g (%g requested)" % (nsamp,mtot,mcluster)

    if mtot > mcluster + tolerance:
        mcum = masses.cumsum()
        last_ind = np.argmax(mcum > mcluster)
        masses = masses[:last_ind]
        mtot = masses.sum()
        if verbose: print "Selected the first %i out of %i masses to get %g total" % (last_ind,len(mcum),mtot)
    else:
        while mtot < mcluster:
            # at least 1 sample, but potentially many more
            nsamp = np.ceil((mcluster-mtot) / mostcommonmass[massfunc])
            newmasses = inverse_imf(np.random.random(nsamp), massfunc=massfunc, **kwargs)
            masses = np.concatenate([masses,newmasses])
            mtot = masses.sum()
            if verbose: print "Sampled %i new stars.  Total is now %g" % (nsamp, mtot)

            if mtot > mcluster+tolerance: # don't force exact equality; that would yield infinite loop
                mcum = masses.cumsum()
                last_ind = np.argmax(mcum > mcluster)
                masses = masses[:last_ind]
                mtot = masses.sum()
                if verbose: print "Selected the first %i out of %i masses to get %g total" % (last_ind,len(mcum),mtot)

    if not silent: print "Total cluster mass is %g (limit was %g)" % (mtot,mcluster)

    return masses

# Vacca Garmany Shull log(lyman continuum) parameters
# Power-law extrapolated from 18 to 8 and from 50 to 150
# (using, e.g., ",".join(["%0.2f" % p for p in polyval(polyfit(log10(vgsmass[:5]),vgslogq[:5],1),log10(linspace(50,150,6)))[::-1]]) 
# where vgsmass does *not* include the extrapolated values)
vgsmass = [150.,  130.,  110.,   90.,   70.,  51.3,44.2,41.0,38.1,35.5,33.1,30.8,28.8,26.9,25.1,23.6,22.1,20.8,19.5,18.4,18.,  16.,  14.,  12.,  10.,   8.][::-1]
vgslogq = [50.51,50.34,50.13,49.88,49.57,49.18,48.99,48.90,48.81,48.72,48.61,48.49,48.34,48.16,47.92,47.63,47.25,46.77,46.23,45.69,45.58,44.65,43.60,42.39,40.96,39.21][::-1]

def lyc_of_star(mass):
    """
    Determine lyman continuum luminosity of a star given its mass
    Uses the Vacca, Garmany, Shull 1996 Table 5 Log Q and Mspec parameters

    returns LogQ
    """

    return np.interp(mass, vgsmass, vgslogq)

def lyc_of_cluster(masses):
    """
    Determine the log of the integrated lyman continuum luminosity of a cluster
    Only M>=8msun count

    masses is a list or array of masses.  
    """
    if max(masses) < 8: return 0
    logq = lyc_of_star(masses[masses >= 8])
    logqtot = np.log10( (10**logq).sum() )
    return logqtot
