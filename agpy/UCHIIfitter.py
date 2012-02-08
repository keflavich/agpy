"""
============
UCHII Fitter
============

Fit a free-free spectrum to an SED.  

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

"""
import pylab
from pylab import *
for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'
try:
    from scipy import optimize
except ImportError:
    print "scipy not installed correctly: UCHIIfitter may fail"
from mpfit import mpfit
import numpy
from agpy.constants import *

#kb = 1.38e-16
#c=3e10
#mu = 1.4
#mh = 1.67e-24
#msun = 1.9889e33
#pc = 3.08568e18      # cm
#au = 1.496e13        # cm
#msun = 1.99e33       # g

unitfactor={'mJy':1e-26,'Jy':1e-23,'cgs':1.0}
freqfactor={'GHz':1e9,'Hz':1.0}


def tnu(Te,nu,EM):
    """
    Te - excitation temperature
    nu - frequency in GHz
    EM - Emission Measure

    Calculates optical depth as a function of temperature, frequency, and
    emission measure from Rohlfs and Wilson 2000's eqns 9.33 and 9.34.

    """
#    nu0 = .3045 * Te**-.643 * EM**.476
    nu0 = Te**1.5 / 1000
    answer_highnu = (nu > nu0) * 3.014e-2 * Te**-1.5 * nu**-2 * EM  
    gff_lownu = ( log(4.955e-2 * nu**-1) + 1.5 * log(Te) )  # <gff> Gaunt factor for free-free
    answer_lownu  = (nu < nu0) * 3.014e-2 * Te**-1.5 * nu**-2 * EM * gff_lownu
    tau = answer_lownu+answer_highnu
    ## altenhoff version
    #tau = 8.235e-2 * Te**-1.35 * nu**-2.1 * EM
    return tau

def Inu(nu,tau,Te,I0=0):
    """
    Calculates flux for a given optical depth, frequency, and temperature
    assuming Rayleigh-Jeans

    nu - frequency in Hz
    tau - optical depth
    Te - excitation temperature (K)
    """
    if I0==0 and isinstance(nu,numpy.ndarray):
        whtau1 = argmin(abs(tau-1))
        nutau1 = nu[whtau1]
        taufactor = 1
    else:
        nutau1 = nu
        taufactor = tau
        """ assumes I0 is set"""
    I0 = 2 * kb * Te * nutau1**2 / c**2 * taufactor
    thin = (tau < 1) * exp(1-tau) * I0
    thick = 2 * kb * Te * (nu * (tau > 1))**2 / c**2
    return thin+thick

def inufit(nu,em,normfac,Te=8500,unit='mJy',frequnit='GHz'):
    """
    Computes the expected intensity as a function of frequency
    for a given emission measure and normalization factor
    nu - array of frequencies (array)
    em - emission measure (float)
    normfac - normalization factor (float)
            - 1/solid angle of source.  1000 AU at 1 kpc = 206265.

    Units: mJy
    """
    _nu = nu*freqfactor[frequnit]
    I0 = 2 * kb * Te * _nu[0]**2 / c**2
    model_intensity = Inu(_nu,tnu(Te,_nu/1e9,em),Te,I0=I0)  # tnu takes GHz
    model_norm = normfac * model_intensity / unitfactor[unit]
    return model_norm


#def inorm(em,nu=freq[1],nu0=freq[0],intens0=flux[0],Te=8500):
#    """
#    Not used?
#    """
#    I0 = 2 * kb * Te * nu0**2 / c**2
#    model_intensity0 = Inu(nu0,tnu(Te,nu0,em),Te,I0=I0)
#    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0)
#    model_norm = intens0/model_intensity0 * model_intensity
#    return model_norm

def inufit_dust(nu,em,normfac,alpha,normfac2,Te=8500):
    """
    inufit with dust added
    """
    I0 = 2 * kb * Te * nu[0]**2 / c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0) 
    model_norm = normfac * model_intensity + normfac2*nu**alpha
    return model_norm

def inufit_dustT(nu,em,normfac,beta,normfac2,dustT,Te=8500):
    I0 = 2 * kb * Te * nu[0]**2 / c**2
    model_intensity = Inu(nu,tnu(Te,nu,em),Te,I0=I0) 
    dustem = 2*hplanck*(nu)**(3+beta) / c**2 * (exp(hplanck*nu*1e9/(kb*abs(dustT))) - 1)**-1
    model_norm = normfac * model_intensity + normfac2/abs(dustT)*dustem
    return model_norm


def mpfitfun(freq,flux,err=None,dust=False,dustT=False):
    """ wrapper around inufit to be passed into mpfit """
    if dust:
        if err == None:
            def f(p,fjac=None): return [0,(flux-inufit_dust(freq,*p))]
        else:
            def f(p,fjac=None): return [0,(flux-inufit_dust(freq,*p))/err]
        return f
    elif dustT:
        if err == None:
            def f(p,fjac=None): return [0,(flux-inufit_dustT(freq,*p))]
        else:
            def f(p,fjac=None): return [0,(flux-inufit_dustT(freq,*p))/err]
        return f
    else:
        if err == None:
            def f(p,fjac=None): return [0,(flux-inufit(freq,*p))]
        else:
            def f(p,fjac=None): return [0,(flux-inufit(freq,*p))/err]
        return f

def emtau(freq,flux,err=None,EMguess=1e7,Te=8500,normfac=5e-6,quiet=1):
    """
    Returns emission measure & optical depth given 
    radio continuum data points at frequency freq with flux
    density flux.

    return bestEM,nu(tau=1),chi^2
    """
    mp = mpfit(mpfitfun(freq,flux,err),xall=[EMguess,normfac],quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror
    chi2 = mp.fnorm
    bestEM = mpp[0]
    normfac = mpp[1]
    nu_tau = ( Te**1.35 / bestEM / 8.235e-2 )**(-1/2.1)

    return bestEM,nu_tau,normfac,chi2

class HIIregion:
    """
    An HII region has properties frequency, flux, and error, which must be
    numpy ndarrays of the same length
    """

    def __init__(self,nu,flux,fluxerr,fluxunit='mJy',frequnit='GHz',beamsize_as2=0.25,dist_kpc=1.0,
            resolved=False,Te=8500,**kwargs):
        order = argsort(asarray(nu))
        self.nu           = asarray(nu)[order]
        self.flux         = asarray(flux)[order]
        self.fluxerr      = asarray(fluxerr)[order]
        self.frequnit     = frequnit
        self.fluxunit     = fluxunit
        self.beamsize_as2 = beamsize_as2
        self.dist_kpc = dist_kpc
        self.resolved = resolved
        self.Te = Te
        self.em,self.nutau,self.normfac,self.chi2 = emtau(self.nu,self.flux,self.fluxerr,Te=self.Te,**kwargs)

    def refit(self,**kwargs):
        """ refit, presumably using different inputs to emtau """
        self.em,self.nutau,self.normfac,self.chi2 = emtau(self.nu,self.flux,self.fluxerr,Te=self.Te,**kwargs)

    def loglogplot(self,numin=1.0,numax=10.0,plottitle='',do_annotations=True,**kwargs):
        x = linspace(numin,numax,500)
        y = inufit(x,self.em,self.normfac)
        loglog(x,y)
        xlabel('Frequency (GHz)')
        ylabel('Flux Density (mJy)')
        title(plottitle)

        errorbar(self.nu,self.flux,yerr=self.fluxerr,fmt=',',**kwargs)

        self.physprops()
        if do_annotations:
            annotate("size (as): %0.2g" % (self.srcsize/au), [.8, .3],textcoords='axes fraction',xycoords='axes fraction')
            annotate("size (au): %0.2g" % (self.srcsize/au), [.8, .3],textcoords='axes fraction',xycoords='axes fraction')
            annotate("mass (msun): %0.2g" % self.mass, [.8, .25],textcoords='axes fraction',xycoords='axes fraction')
            annotate("EM: %0.2g" % self.em, [.8, .2],textcoords='axes fraction',xycoords='axes fraction')
            annotate("Nu(Tau=1): %0.2g" % self.nutau, [.8, .15],textcoords='axes fraction',xycoords='axes fraction')
            annotate("N(lyc): %0.2g" % self.Nlyc, [.8,.1],textcoords='axes fraction',xycoords='axes fraction')
            annotate("dens: %0.2g" % self.dens, [.8,.05],textcoords='axes fraction',xycoords='axes fraction')

    def physprops(self):
        """
        Get the source size (au), density (cm^-3), 
        mass (msun), and Nlyc of the UCHII

        Also return EM and nutau

        ERROR IN CURRENT VERSION
        """
        if self.resolved:
            self.srcsize = self.beamsize_as2 * (self.dist_kpc*1000.0*au)**2
        else:
            self.srcsize = sqrt(self.flux[0]*unitfactor[self.fluxunit]/(2*kb*self.Te) * \
                    (c/(self.nu[0]*freqfactor[self.frequnit]))**2 * (self.dist_kpc*1e3*pc)**2 / pi) 
        self.dens = sqrt(self.em/(self.srcsize/pc))
        self.mass = self.dens * 4.0/3.0 * pi * self.srcsize**3 * muh * mh / msun

        U = self.dens**(2/3.) * self.srcsize/pc
        self.Nlyc = 8.04e46*self.Te**-.85 * U**3

        return self.srcsize/au,self.dens,self.mass,self.Nlyc,self.em,self.nutau



    # Cara test data:
    # nu = array([1.4,5,8.33]); flux=array([4.7,9.2,9.1]); err=array([.52,.24,.07])
    # em,nutau,normfac,chi2 = UCHIIfitter.emtau(nu,flux,err)

__all__ = [tnu,Inu,unitfactor,freqfactor,inufit,emtau,mpfitfun,HIIregion]
