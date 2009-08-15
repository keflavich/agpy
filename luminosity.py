
import numpy
from numpy import pi
pc = 3.086e18 # cm
c = 2.99792458e10 # cm/s
lsun = 3.839e33 # erg/s


class luminosity:
    """
    Measures luminosity from an SED following the method of Kauffmann et al 2008

    For frequency, bandwidth, etc. see http://casa.colorado.edu/~ginsbura/filtersets.htm

    Units are CGS

    nu - frequency (assumed Hz)
    wnu - width of frequency bin
    fnu - flux (assumed Jy)
    efnu - flux error
    """

    def __init__(self,nu,fnu,wnu=None,efnu=None,dist_pc=None,npoints=1e5):

        self.nu   = numpy.asarray(nu)
        nusort = numpy.argsort(self.nu)
        self.nu = self.nu[nusort]
        try:
            self.wnu  = numpy.asarray(wnu)[nusort]
        except:
            self.wnu = wnu
        self.fnu  = numpy.asarray(fnu)[nusort]
        try: 
            self.efnu = numpy.asarray(efnu)[nusort]
        except:
            self.efnu = efnu
        self.dist_pc = dist_pc
        self.init_interp(npoints)

    def fbol_meas(self):
        """
        Returns the total integrated flux (int[nu fnu dnu]) measured within 
        the bands.  Does not interpolate

        The Casoli et al 1986 formula is not used...
        FIR (10^-13 W m^-2) = 1.75 ( F12 / 0.79 + F25/2 + F60/3.9 + F100/9.9 )
        """
        self._intflux = ( self.fnu * self.wnu ).sum()
        return self._intflux

    def lbol_meas(self,dist_pc=None):
        """
        Returned the luminosity within the measured bins.  Does not interpolate.
        A distance must be specified either here or in the initialization
        """

        if dist_pc is not None:
            self.dist_pc = dist_pc
        elif self.dist_pc is None:
            raise ValueError("Must specify distance to compute luminosity")

        lum = 4*pi*(self.dist_pc*pc)**2 * self.fbol_meas() * 1e-23 / lsun

        return lum

    def init_interp(self,npoints):
        """
        initializes the interpolation parameters
        """
        self.interpnu = 10**numpy.linspace(numpy.log10(self.nu[0]),numpy.log10(self.nu[-1]),npoints)
        self.interpfnu = self.interpnu*0.0
        self.addedpoint = False

    def mminterp(self,freq,lowfreq=c/500e-4,alpha=4.0,addpoint=True):
        """
        Creates an interpolated point using an assumed nu^alpha opacity...
        defaults to alpha=4

        Default point location is 500 microns
        """

        if freq is None:
            whpt = numpy.argmin(self.nu)
        else:
            whpt = (self.nu==freq)

        if whpt == numpy.array([]):
            raise ValueError("Couldn't find a data point at that frequency.  Please check your \
                    input frequency, input data, or a floating point match.")

        numm = self.nu[whpt]
        fnumm = self.fnu[whpt]

        userange = (self.interpnu < freq) * (self.interpnu > lowfreq)
        self.interpfnu[userange] = fnumm * (self.interpnu[userange]/numm)**alpha

        if addpoint and not self.addedpoint:
            newfnu = fnumm * (lowfreq/numm)**alpha
            self.fnu = numpy.concatenate((self.fnu,numpy.array([newfnu])))
            self.nu = numpy.concatenate((self.nu,numpy.array([lowfreq])))
            nusort = numpy.argsort(self.nu)
            self.nu = self.nu[nusort]
            self.fnu = self.fnu[nusort]
            self.addedpoint = True


    def fbol_interp(self,mminterp=True,npoints=1e5,addpoint=True,mmfreq=None,extrap=False):
        """
        Interpolates between data points to integrate over SED

        If mminterp is set, will assume a nu^4 power law (opacity lambda^-2)

        Returns int( nuFnu ) in units HzJy

        Extrapolate via a constant line at the low/high frequency data point
        """

        if mminterp:
            self.mminterp(mmfreq,addpoint=addpoint)

        # powerlaw interpolation is done in logspace
        logf = numpy.log10(self.fnu)
        logn = numpy.log10(self.nu)

        # slope = dy/dx
        alpha = ( logf[:-1]-logf[1:] ) / ( logn[:-1]-logn[1:] )

        # evenly spaced sampling in logspace
        newnu  = self.interpnu
        newfnu = self.interpfnu

        # need to loop through data and interpolate
        for ind in numpy.arange(len(logn)-1):
            lownu  = self.nu[ind]
            highnu = self.nu[ind+1]
            userange = (newnu > lownu) * (newnu < highnu)
            nu0 = newnu[userange][0]
            newfnu[userange] = self.fnu[ind] * (newnu[userange]/nu0)**alpha[ind]

        dnu = newnu*0
        dnu[:-1] = newnu[:-1]-newnu[1:]
        dnu[-1] = dnu[-2]
        self.interpdnu = numpy.abs(dnu)

        self.interpnu = newnu
        self.interpfnu = newfnu

        self._fbol_integ = (self.interpdnu*self.interpfnu).sum()

        return self._fbol_integ

    def lbol_interp(self,**kwargs):
        """
        Bolometric luminosity from interpolation in units of solar luminosities
        """

        self._lbol_interp = 4*pi*(self.dist_pc*pc)**2 * self.fbol_interp(**kwargs) * 1e-23 / lsun

        return self._lbol_interp

    """ 
    Test case:

    nu = c/array([12e-4,25e-4,60e-4,100e-4,.12])
    fnu = [4.97,26.08,66.34,80.04,3.06]
    """

import readcol

def test_case():
    """ A specific test case using the Klein 2005 luminosity computations """

    dtab = readcol.readcol('/Users/adam/agpy/tests/klein2005sourcelist.txt',skipline=33,fsep='|')

    kleinlum = dtab[:,5].astype('float')

    irasfluxes = dtab[:,numpy.array([22,25,28,31])].astype('float')
    dist = dtab[:,3].astype('float')*1000

    mylum = kleinlum*0
    mylum_interp = kleinlum*0
    nu = c*1e4/numpy.array([12.0,25.0,60.0,100.0])
    dnu = numpy.array([1.53e13,5.789e12,3.75e12, 1.114e12])

    mlarr = []

    for ind in xrange(len(mylum)):
        ml = luminosity(nu,irasfluxes[ind],dnu,dist_pc=dist[ind])
        mylum[ind] = ml.lbol_meas()
        mylum_interp[ind] = ml.lbol_interp(addpoint=False)
        mlarr.append(ml)

    return kleinlum,mylum,mylum_interp,mlarr







