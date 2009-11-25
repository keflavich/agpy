
import numpy
try: 
    import pylab
    pylabok = True
except:
    pylabok = False
    print "Pylab could not be imported; plotting functions will be disabled"
try:
    import scipy.special
    scipyok = True
except:
    scipyok = False
    print "Scipy could not be imported."
from numpy import pi
pc = 3.086e18 # cm
c = 2.99792458e10 # cm/s
lsun = 3.839e33 # erg/s
kb = 1.3806503e-16 # erg/K
h = 6.626068e-27 # erg s


class luminosity:
    """
    Measures luminosity from an SED following the method of Kauffmann et al 2008
    http://adsabs.harvard.edu/abs/2008A%26A...487..993K

    For frequency, bandwidth, etc. see http://casa.colorado.edu/~ginsbura/filtersets.htm

    Units are CGS

    nu - frequency (assumed Hz)
    wnu - [optional/required for non-interpolation] width of frequency bin
    lnu/unu - lower/upper bounds on frequency bin
    fnu - flux (assumed Jy)
    efnu - [optional; currently does nothing] flux error

    When interpolation is used, data is added to the luminosity' classes nu/fnu vectors.
    If you use, e.g. luminosity.lbol_meas(), it will still work correctly because it sets
    wnu ~ unu-lnu = 0
    """

    def __init__(self,nu,fnu,wnu=None,lnu=None,unu=None,efnu=None,dist_pc=None,npoints=1e5):

        self.nu   = numpy.asarray(nu)
        nusort = numpy.argsort(self.nu)
        self.nu = self.nu[nusort]
        try:
            self.wnu  = numpy.asarray(wnu)[nusort]
        except:
            self.wnu = wnu
        try:
            self.lnu  = numpy.asarray(lnu)[nusort]
            self.unu  = numpy.asarray(unu)[nusort]
        except:
            self.lnu = lnu
            self.unu = unu
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
        if self.wnu is None and self.lnu is None:
            print "Must specify bandwidths."
            return
        elif self.wnu is not None:
            self._intflux = ( self.fnu * self.wnu ).sum()
        elif self.lnu is not None and self.unu is not None:
            self._intflux = ( self.fnu * (self.unu-self.lnu) ).sum()
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
        self.extrapolated = False

    def mminterp(self,freq,lowfreq=c/500e-4,alpha=4.0,addpoint=True):
        """
        Creates an interpolated point using an assumed nu^alpha opacity...
        defaults to alpha=4 (beta=2)

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

            if self.efnu is not None:
                self.efnu = numpy.concatenate((self.efnu,numpy.array([0.0])))
                self.efnu = self.efnu[nusort]
            if self.wnu is not None:
                self.wnu = numpy.concatenate((self.wnu,numpy.array([0.0])))
                self.wnu = self.wnu[nusort]
            if self.lnu is not None:
                self.lnu = numpy.concatenate((self.lnu,numpy.array([lowfreq])))
                self.lnu = self.lnu[nusort]
            if self.unu is not None:
                self.unu = numpy.concatenate((self.unu,numpy.array([lowfreq])))
                self.unu = self.unu[nusort]


    def fbol_interp(self,fnu=None,mminterp=True,npoints=1e5,addpoint=True,mmfreq=None,extrap=True,write=True):
        """
        Interpolates between data points to integrate over SED

        If mminterp is set, will assume a nu^4 power law (opacity lambda^-2)
        from the longest wavelength point

        Returns int( nuFnu ) in units HzJy

        Extrapolate via a constant line at the low/high frequency data point
        """

        if mminterp:
            self.mminterp(mmfreq,addpoint=addpoint)

        if fnu is None:
            fnu = self.fnu

        # powerlaw interpolation is done in logspace
        logf = numpy.log10(fnu)
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
            newfnu[userange] = fnu[ind] * (newnu[userange]/nu0)**alpha[ind]

        dnu = newnu*0
        dnu[:-1] = newnu[:-1]-newnu[1:]
        dnu[-1] = dnu[-2]

        if extrap and not self.extrapolated:
            # extrapolate 50% from either endpoint
            nulow = newnu.min()/2.0
            nuhigh = newnu.max()*2.0
            dnulow = newnu.min()-nulow
            dnuhigh = newnu.max()-nuhigh
            scale_lownu  = min( [ newfnu[0]/newfnu[1]  , 1.0 ] ) # do not allow positive-slope extrapolation
            scale_highnu = min( [ self.fnu[-2]/self.fnu[-1], 0.5 ] )
            newnu  = numpy.concatenate( ([nulow],newnu,[nuhigh]) )
            newfnu = numpy.concatenate( ([newfnu[0]*scale_lownu],newfnu,[newfnu[-1]*scale_highnu]) )
            dnu    = numpy.concatenate( ([dnulow],dnu,[dnuhigh]) )
            self.extrapolated=True

        if write:
            self.interpdnu = numpy.abs(dnu)
            self.interpnu = newnu
            self.interpfnu = newfnu

            fbint = (self.interpdnu*self.interpfnu).sum()

            self._fbol_integ = fbint
        else: # unfortunate repeated code; should clean this up / not write self. variables
            interpdnu = numpy.abs(dnu)
            interpnu = newnu
            interpfnu = newfnu

            fbint = (interpdnu*interpfnu).sum()

        return fbint

    def lbol_interp(self,**kwargs):
        """
        Bolometric luminosity from interpolation in units of solar luminosities

        By default, adds a point at 500 microns extrapolated using a nu^4 power law 
        (opacity nu^2) from the longest wavelength data point.  Specify addpoint=False
        to disable this feature.
        """

        self._lbol_interp = 4*pi*(self.dist_pc*pc)**2 * self.fbol_interp(**kwargs) * 1e-23 / lsun

        return self._lbol_interp

    def lbol_interp_ulim(self,**kwargs):
        """
        If errors are specified, returns lbol(fnu+efnu)
        """
        if self.efnu is None:
            print "Can't calculate limits unless errors are specified"

        lbol_upper = 4*pi*(self.dist_pc*pc)**2 * self.fbol_interp(fnu=self.fnu+self.efnu,write=False,**kwargs) * 1e-23 / lsun

        return lbol_upper

    def lbol_interp_llim(self,**kwargs):
        """
        If errors are specified, returns lbol(fnu+efnu)
        """
        if self.efnu is None:
            print "Can't calculate limits unless errors are specified"

        lbol_upper = 4*pi*(self.dist_pc*pc)**2 * self.fbol_interp(fnu=self.fnu-self.efnu,write=False,**kwargs) * 1e-23 / lsun

        return lbol_upper

    def tbol(self,interp=True):
        """
        Computes the "bolometric temperature" as specified in the same document.
        Uses scipy's zeta function if scipy is available
        
        """
        #print "tbol not yet implemented"

        if scipyok:
            zeta4d5 = scipy.special.zeta(4,1)/scipy.special.zeta(5,1)
        else:
            zeta4d5 = 1.0437788248434832

        if self.interpdnu is not None and interp:
            meannu = (self.interpnu*self.interpfnu*self.interpdnu).sum() / (self.interpfnu*self.interpdnu).sum()
        elif self.lnu is not None and self.unu is not None:
            meannu = (self.nu*self.fnu*(self.unu-self.lnu)).sum() / (self.fnu*(self.unu-self.lnu)).sum()
        elif self.wnu is not None:
            meannu = (self.nu*self.fnu*self.wnu).sum() / (self.fnu*self.wnu).sum()
        else:
            print "Need to specify either wnu (width of each band) or interpolate."
            return

        self.Tbol = zeta4d5 * h * meannu / (4 * kb)

        return self.Tbol


    def plotsed(self,loglog=True,nufnu=False,interpplot=False,**kwargs):
        """ Plots the SED """
        if not pylabok:
            print "pylab was not successfully imported.  Aborting."
            return

        # must divide wnu by 2 because plotter assumes "1-sigma", but wnu is full width

        if self.efnu is None:
            efnu = 0
        else:
            efnu = self.efnu

        if nufnu:
          if self.lnu is not None and self.unu is not None:
            xerr = numpy.abs( numpy.array([self.lnu,self.unu]) - self.nu )
            pylab.errorbar(self.nu,self.nu*self.fnu,xerr=xerr,yerr=efnu*self.nu,fmt=',',**kwargs)
          elif self.wnu is not None:
            pylab.errorbar(self.nu,self.nu*self.fnu,xerr=self.wnu/2.0,yerr=efnu*self.nu,fmt=',',**kwargs)
          else:
            pylab.errorbar(self.nu,self.nu*self.fnu,yerr=efnu*self.nu,fmt=',',**kwargs)
          if interpplot:
            pylab.plot(self.interpnu,self.interpfnu*self.interpnu,**kwargs)
        else:
          if self.lnu is not None and self.unu is not None:
            xerr = numpy.abs( numpy.array([self.lnu,self.unu]) - self.nu )
            pylab.errorbar(self.nu,self.fnu,xerr=xerr,yerr=efnu,fmt=',',**kwargs)
          elif self.wnu is not None:
            pylab.errorbar(self.nu,self.fnu,xerr=self.wnu/2.0,yerr=efnu,fmt=',',**kwargs)
          else:
            pylab.errorbar(self.nu,self.fnu,yerr=efnu,fmt=',',**kwargs)
          if interpplot:
            pylab.plot(self.interpnu,self.interpfnu,**kwargs)
        ax = pylab.gca()
        if loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
        else:
            ax.set_xscale('linear')
            ax.set_yscale('linear')
        pylab.xlabel('Frequency (Hz)')
        pylab.ylabel('Flux Density (Jy)')

        return ax


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
        mylum_interp[ind] = ml.lbol_interp(addpoint=False,extrap=True)
        mlarr.append(ml)

    return kleinlum,mylum,mylum_interp,mlarr

""" 
Test case:

nu = c/array([12e-4,25e-4,60e-4,100e-4,.12])
fnu = [4.97,26.08,66.34,80.04,3.06]

# G31.28
nu = c/array([8.3e-4,12e-4,21.3e-4,25e-4,60e-4,100e-4,.045,.085,.12])
fnu = array([2.4,4.85,34.9,89.33,1071,3693,160,19.5,5.4])

"""

