"""
============================
Simple black-body calculator
============================

Includes both wavelength and frequency blackbody functions.  Has flexible
units.  Also allows for a few varieties of modified blackbody.

"""
try:
    from numpy import exp
except ImportError:
    from math import exp

unitdict = {'cgs':{'h':6.626068e-27,
                   'k':1.3806503e-16,
                   'c':2.99792458e10,
                   'mh':1.67262158e-24 * 1.00794,
                   'length':'cm'},
            'mks':{'h':6.626068e-34,
                   'k':1.3806503e-23,
                   'c':2.99792458e8,
                   'mh':1.67262158e-27 * 1.00794,
                   'length':'m'}
            }

frequency_dict = {'Hz':1.0,
                  'kHz':1e3,
                  'MHz':1e6,
                  'GHz':1e9,
                  'THz':1e12,
                  }

def blackbody(nu,temperature, scale=1.0, units='cgs',frequency_units='Hz',
              normalize=max, beta=0):
    # load constants in desired units
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']

    # convert nu to Hz
    nu = nu * frequency_dict[frequency_units]

    I = 2*h*nu**3 / c**2 * (exp(h*nu/(k*temperature)) - 1)**-1

    if normalize and hasattr(I,'__len__'):
        if len(I) > 1:
            return I/normalize(I) * scale
        else:
            return I * scale
    else:
        return I * scale

wavelength_dict = {'meters':1.0,'m':1.0,
                   'centimeters':1e-2,'cm':1e-2,
                   'millimeters':1e-3,'mm':1e-3,
                   'nanometers':1e-9,'nm':1e-9,
                   'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
                   'kilometers':1e3,'km':1e3,
                   'angstroms':1e-10,'A':1e-10,'Angstroms':1e-10,
                   }

def blackbody_wavelength(lam,temperature, scale=1.0,
        units='cgs',wavelength_units='Angstroms', normalize=max, beta=0):
    # load constants in desired units
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']

    # converta lambd to cm/m
    lam = lam * wavelength_dict[wavelength_units] / (1e-2 if units=='cgs' else 1)

    I = 2*h*c**2 / lam**5 * (exp(h*c/(k*temperature*lam)) - 1)**-1

    if normalize and hasattr(I,'__len__'):
        if len(I) > 1:
            return I/normalize(I) * scale
        else:
            return I * scale
    else:
        return I * scale

def modified_blackbody(nu, temperature, beta=1.75, logN=22, logscale=0.0,
                       muh2=2.8, units='cgs',frequency_units='Hz', kappa0=4.0,
                       nu0=505e9, normalize=max, dusttogas=100.):
    """
    Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1  (1 - e^(-tau_nu) )
    Kappa0 and Nu0 are set as per http://arxiv.org/abs/1101.4654 which uses OH94 values.
    beta = 1.75 is a reasonable default for Herschel data
    N = 1e22 is the column density in cm^-2

    nu0 and nu must have same units!

    Parameters
    ----------
    nu : float
        Frequency in units of `frequency_units`
    temperature : float
        Temperature in Kelvins
    beta : float
        The blackbody modification value; the blackbody function is multiplied
        by :math:`(1-exp(-(\\nu/\\nu_0)**\\beta))`
    logN : float
        The log column density to be fit
    logscale : float
        An arbitrary logarithmic scale to apply to the blackbody function
        before passing it to mpfit; this is meant to prevent numerical
        instability when attempting to fit very small numbers.
        Can also be used to represent, e.g., steradians
    muh2 : float
        The mass (in amu) per molecule of H2.  Defaults to 2.8.
    units : 'cgs' or 'mks'
        The unit system to use
    frequency_units : string
        Hz or some variant (GHz, kHz, etc)
    kappa0 : float
        The opacity in cm^2/g *for gas* at nu0 (see dusttogas)
    nu0 : float
        The frequency at which the opacity power law is locked
    normalize : function or None
        A normalization function for the blackbody.  Set to None if you're
        interested in the amplitude of the blackbody
    dusttogas : float
        The dust to gas ratio.  The opacity kappa0 is divided by this number to
        get the opacity of the dust
    """
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']
    mh = unitdict[units]['mh']

    kappanu = kappa0 / dusttogas * (nu/nu0)**beta
    # numpy apparently can't multiply floats and longs
    tau = muh2 * mh * kappanu * 10.0**logN

    modification = (1.0 - exp(-1.0 * tau))

    I = blackbody(nu, temperature, units=units,
                  frequency_units=frequency_units,
                  normalize=normalize)*modification

    if normalize and hasattr(I,'__len__'):
        if len(I) > 1:
            return I/normalize(I) * 10.**logscale
        else:
            return I * 10.**logscale
    else:
        return I * 10.**logscale

def greybody(nu, temperature, beta, A=1.0, logscale=0.0,
             units='cgs', frequency_units='Hz', kappa0=4.0, nu0=3000e9,
             normalize=max):
    """
    Same as modified blackbody... not sure why I have it at all, though the
    normalization constants are different.
    """
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']

    modification = (1. - exp(-(nu/nu0)**beta))
    I = blackbody(nu,temperature,units=units,frequency_units=frequency_units,normalize=normalize)*modification

    if normalize and hasattr(I,'__len__'):
        if len(I) > 1:
            return I/normalize(I) * 10.**logscale
        else:
            return I * 10.**logscale
    else:
        return I * 10.**logscale

def modified_blackbody_wavelength(lam, temperature, beta=1.75, logN=22,
                                  logscale=0.0, muh2=2.8, units='cgs',
                                  wavelength_units='Angstroms', kappa0=4.0,
                                  nu0=505e9, dusttogas=100., normalize=max):
    """
    Snu =  2hnu^3 c^-2  (e^(hnu/kT) - 1)^-1  (1 - e^(-tau_nu) )
    Kappa0 and Nu0 are set as per http://arxiv.org/abs/1101.4654 which uses OH94 values.
    beta = 1.75 is a reasonable default for Herschel data
    N = 1e22 is the column density in cm^-2

    nu0 and nu must have same units!  But wavelength is converted to frequency
    of the right unit anyway

    Parameters
    ----------
    lam : float
        Wavelength in units of `wavelength_units`
    temperature : float
        Temperature in Kelvins
    beta : float
        The blackbody modification value; the blackbody function is multiplied
        by :math:`(1-exp(-(\\nu/\\nu_0)**\\beta))`
    logN : float
        The log column denisty to be fit
    logscale : float
        An arbitrary logarithmic scale to apply to the blackbody function
        before passing it to mpfit; this is meant to prevent numerical
        instability when attempting to fit very small numbers.
        Can also be used to represent, e.g., steradians
    muh2 : float
        The mass (in amu) per molecule of H2.  Defaults to 2.8.
    units : 'cgs' or 'mks'
        The unit system to use
    wavelength_units : string
        A valid wavelength (e.g., 'angstroms', 'cm','m')
    kappa0 : float
        The opacity in cm^2/g *for gas* at nu0 (see dusttogas)
    nu0 : float
        The frequency at which the opacity power law is locked.
        kappa(nu) = kappa0/dusttogas * (nu/nu0)**beta
    normalize : function or None
        A normalization function for the blackbody.  Set to None if you're
        interested in the amplitude of the blackbody
    dusttogas : float
        The dust to gas ratio.  The opacity kappa0 is divided by this number to
        get the opacity of the dust
    """
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']
    mh = unitdict[units]['mh']

    nu = c/(lam*wavelength_dict[wavelength_units]/wavelength_dict[unitdict[units]['length']])

    #I = modified_blackbody(nu, temperature, beta=beta, frequency_units='Hz',
    #                       normalize=normalize, nu0=nu0, kappa0=kappa0,
    #                       muh2=muh2, logscale=logscale, units=units,
    #                       logN=logN)

    kappanu = kappa0/dusttogas * (nu/nu0)**beta
    tau = muh2 * mh * kappanu * 10.**logN

    modification = (1.0 - exp(-1.0 * tau))

    I = blackbody(nu,temperature,units=units,frequency_units='Hz',normalize=normalize)*modification

    if normalize and hasattr(I,'__len__'):
        if len(I) > 1:
            return I/normalize(I) * 10.**logscale
        else:
            return I * 10.**logscale
    else:
        return I * 10**logscale


try:
    import agpy.mpfit as mpfit

    def fit_blackbody(xdata, flux, guesses=(0,0), err=None,
            blackbody_function=blackbody, quiet=True, **kwargs):
        """
        Parameters
        ----------
        xdata : array
            Array of the X-values (frequency, wavelength) of the data
        flux : array
            The fluxes corresponding to the xdata values
        guesses : (Temperature,Scale) or (Temperature,Beta,Scale)
            The input guesses.  3 parameters are used for greybody 
            fitting, two for temperature fitting.
        blackbody_function: function
            Must take x-axis (e.g. frequency), temperature, scale, and then
            optionally beta args 
        quiet : bool
            quiet flag passed to mpfit

        Returns
        -------
        mp : mpfit structure
            An mpfit structure.  Access parameters and errors via
            `mp.params` and `mp.perror`.  The covariance matrix
            is in mp.covar.

        Examples
        --------
        >>> wavelength = array([20,70,160,250,350,500,850,1100])
        >>> flux = modified_blackbody_wavelength(wavelength, 15, beta=1.75,
                logN=22, wavelength_units='microns', normalize=False,
                logscale=16)
        >>> err = 0.1 * flux
        >>> np.random.seed(0)
        >>> flux += np.random.randn(len(wavelength)) * err
        >>> tguess, bguess, nguess = 20.,2.,21.5
        >>> mp = fit_blackbody(wavelength, flux, err=err,
                 blackbody_function=modified_blackbody_wavelength,
                 logscale=16, guesses=(tguess, bguess, nguess),
                 wavelength_units='microns')
        >>> print mp.params 
        [ 14.99095224   1.78620237  22.05271119]
        >>> # T~14.9 K, beta ~1.79, column ~10^22
        """

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-blackbody_function(x, *p,
                    normalize=False, **kwargs))]
            else:
                def f(p,fjac=None): return [0,(y-blackbody_function(x, *p,
                    normalize=False, **kwargs))/err]
            return f

        err = err if err is not None else flux*0.0 + 1.0

        mp = mpfit.mpfit(mpfitfun(xdata,flux,err), guesses, quiet=quiet)

        return mp
except ImportError:
    pass

try:
    import lmfit
    try:
        from collections import OrderedDict
    except ImportError:
        from ordereddict import OrderedDict
    import numpy as np
    def fit_blackbody_lmfit(xdata, flux, guesses=(0,0), err=None,
            blackbody_function=blackbody, quiet=True, **kwargs):
        """
        Parameters
        ----------
        xdata : array
            Array of the X-values (frequency, wavelength) of the data
        flux : array
            The fluxes corresponding to the xdata values
        guesses : (Temperature,Scale) or (Temperature,Beta,Scale)
            The input guesses.  3 parameters are used for greybody 
            fitting, two for temperature fitting.
        blackbody_function: function
            Must take x-axis (e.g. frequency), temperature, scale, and then
            optionally beta args 
        quiet : bool
            quiet flag passed to mpfit

        kwargs are past to blackbody function

        Examples
        --------
        >>> wavelength = np.array([20,70,160,250,350,500,850,1100])
        >>> flux = modified_blackbody_wavelength(wavelength, 15, beta=1.75,
                wavelength_units='microns', normalize=False, logN=22, logscale=16)
        >>> err = 0.1 * flux
        >>> flux += np.random.randn(len(wavelength)) * err
        >>> tguess, bguess, nguess = 20.,2.,21.5
        >>> lm = fit_blackbody_lmfit(wavelength, flux, err=err,
                 blackbody_function=modified_blackbody_wavelength, logscale=16,
                 guesses=(tguess,bguess,nguess),
                 wavelength_units='microns')
        >>> print lm.params
        
        >>> # If you want to fit for a fixed beta, do this:
        >>> parameters = lmfit.Parameters(OrderedDict([ (n,lmfit.Parameter(x)) for n,x
                in zip(('T','beta','N'),(20.,2.,21.5)) ]))
        >>> import lmfit
        >>> parameters['beta'].vary = False
        >>> lm = fit_blackbody_lmfit(wavelength, flux, err=err,
                 blackbody_function=modified_blackbody_wavelength, logscale=16,
                 guesses=parameters,
                 wavelength_units='microns')
        >>> print lm.params
        """

        def lmfitfun(x,y,err):
            if err is None:
                def f(p): return (y-blackbody_function(x, *[p[par].value for par in p],
                    normalize=False, **kwargs))
            else:
                def f(p): return (y-blackbody_function(x, *[p[par].value for par in p],
                    normalize=False, **kwargs))/err
            return f

        if not isinstance(guesses,lmfit.Parameters):
            guesspars = lmfit.Parameters(
                    OrderedDict([ (n,lmfit.Parameter(value=x,name=n))
                        for n,x in zip(('T','beta','N'),guesses) ]))
        else:
            guesspars = guesses

        minimizer = lmfit.minimize( lmfitfun(xdata,np.array(flux),err),
                guesspars)

        return minimizer
except ImportError:
    pass

try:
    import pymodelfit

    # FAILS:
    # SyntaxError: can't have kwargs in model function
    #class pmf_blackbody(pymodelfit.FunctionModel1DAuto):
    #    def f(self, x, T=20.0, scale=1.0, beta=1.5,
    #            blackbody_function=blackbody, **kwargs):
    #        return blackbody_function(x, T, scale, beta=beta)
except ImportError:
    pass

try:
    import numpy as np
    old_errsettings = np.geterr()
    import pymc # pymc breaks np error settings
    np.seterr(**old_errsettings)

    def fit_blackbody_montecarlo(frequency, flux, err=None,
            temperature_guess=10, beta_guess=None, scale_guess=None,
            blackbody_function=blackbody, quiet=True, return_MC=False,
            nsamples=5000, burn=1000, min_temperature=0, max_temperature=100,
            scale_keyword='scale', max_scale=1e60,
            multivariate=False,
            **kwargs):
        """
        Parameters
        ----------
        frequency : array
            Array of frequency values
        flux : array
            array of flux values
        err : array (optional)
            Array of error values (1-sigma, normal)
        temperature_guess : float
            Input / starting point for temperature
        min_temperature : float
        max_temperature : float
            Lower/Upper limits on fitted temperature
        beta_guess : float (optional)
            Opacity beta value
        scale_guess : float
            Arbitrary scale value to apply to model to get correct answer
        blackbody_function: function
            Must take x-axis (e.g. frequency), temperature, then scale and beta
            keywords (dependence on beta can be none)
        return_MC : bool
            Return the pymc.MCMC object?
        nsamples : int
            Number of samples to use in determining the posterior distribution
            (the answer)
        burn : int
            number of initial samples to ignore
        scale_keyword : ['scale','logscale','logN']
            What scale keyword to pass to the blackbody function to determine
            the amplitude
        kwargs : kwargs
            passed to blackbody function
        """

        d = {}
        d['temperature'] = pymc.distributions.Uniform('temperature',
                min_temperature, max_temperature, value=temperature_guess)
        d['scale'] = pymc.distributions.Uniform('scale',0,max_scale,
                value=scale_guess)
        if beta_guess is not None:
            d['beta'] = pymc.distributions.Uniform('beta',0,10,
                    value=beta_guess)
        else:
            d['beta'] = pymc.distributions.Uniform('beta',0,0,
                    value=0)


        @pymc.deterministic
        def bb_model(temperature=d['temperature'],
                scale=d['scale'],
                beta=d['beta']):
            kwargs[scale_keyword] = scale
            y = blackbody_function(frequency, temperature, 
                    beta=beta, normalize=False, **kwargs)
            #print kwargs,beta,temperature,(-((y-flux)**2)).sum()
            return y

        d['bb_model'] = bb_model

        if err is None:
            d['err'] = pymc.distributions.Uninformative('error',value=1.)
        else:
            d['err'] = pymc.distributions.Uninformative('error',value=err,observed=True)

        d['flux'] = pymc.distributions.Normal('flux', mu=d['bb_model'], tau=1./d['err']**2,
                value=flux, observed=True)

        #print d.keys()
        MC = pymc.MCMC(d)
        
        if nsamples > 0:
            MC.sample(nsamples, burn=burn)
            if return_MC:
                return MC

            MCfit = pymc.MAP(MC)
            MCfit.fit()
            T = MCfit.temperature.value
            scale = MCfit.scale.value

            if beta_guess is not None:
                beta = MCfit.beta.value
                return T,scale,beta
            else:
                return T,scale

        return MC

except ImportError:
    pass


if __name__=="__main__":

    print "Fitting tests"
    import itertools
    import numpy as np
    import agpy
    import pylab
    import matplotlib
    temperatures = [5,10,15,20,25]
    betas        = [1,1.25,1.5,1.75,2.0,2.25,2.5]
    columns      = [22,23,24]
    wavelengths = [#np.array([70,160,250,350,500,850,1100],dtype='float'),
        np.array([160,250,350,500,1100],dtype='float'),
        np.array([160,250,350,500],dtype='float')]
    errlevels = [0.1,0.2,0.05]

    temperature = 20.
    beta = 1.75
    column = 22
    wavelength = wavelengths[0]
    errlevel = errlevels[0]

    #for temperature,beta,column in itertools.product(temperatures,betas,columns):

    tguess=20
    bguess=1.5
    nguess=21
    tguess=temperature
    bguess=beta
    nguess=column

    bbmcs = {}

    MCtest = False
    
    if MCtest: 
        for ii,(errlevel,wavelength) in enumerate(itertools.product(errlevels,wavelengths)):
            flux = modified_blackbody_wavelength(wavelength, temperature,
                    beta=beta, wavelength_units='microns', normalize=False, logN=column,
                    logscale=16)
            err = flux*errlevel

            bbmc = fit_blackbody_montecarlo(wavelength, flux,
                    blackbody_function=modified_blackbody_wavelength,
                    return_MC=True, wavelength_units='microns', nsamples=1,
                    scale_guess=nguess, beta_guess=bguess, temperature_guess=tguess,
                    scale_keyword='logN', max_scale=30, err=err, logscale=16,
                    burn=0)
            flux = bbmc.flux.rand()
            bbmc = fit_blackbody_montecarlo(wavelength, flux,
                    blackbody_function=modified_blackbody_wavelength,
                    return_MC=True, wavelength_units='microns', nsamples=30000,
                    scale_guess=nguess, beta_guess=bguess, temperature_guess=tguess,
                    scale_keyword='logN', max_scale=30, err=err, logscale=16,
                    burn=10000)
            bbmcs[errlevel+len(wavelength)] = bbmc
            mp = fit_blackbody(wavelength, flux, err=err,
                    blackbody_function=modified_blackbody_wavelength,
                    logscale=16, guesses=(tguess, bguess, nguess), wavelength_units='microns')
            print
            print "      %10s %10s %10s %10s %10s %10s %10s %10s %10s S/N=%f wls: %s" % ('T','e(T)','B','e(B)','N','e(N)','T-B','B-N','T-N',1/errlevel,wavelength)
            print "input %10.3g %10s %10.3g %10s %10.3g %10s " % (temperature,"",beta,"",column,"")
            mppars = [p for a in zip(mp.params,mp.perror) for p in a]
            print("chi^2 %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g " % tuple(mppars)+
                  "%10.3g %10.3g %10.3g " % (mp.covar[1,0],mp.covar[2,1],mp.covar[2,0]))
            mcstats = bbmc.stats()
            try:
                N = pymc.NormApprox(bbmc)
                N.fit()
                Ncov = N.C[N.temperature,N.beta,N.scale]
                mcpars = [p for a in zip(N.mu[N.temperature, N.beta, N.scale],
                                         np.array(Ncov.diagonal())[0]) for p in a]
                print("MCMC  %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g " % tuple(mcpars)+
                      "%10.3g %10.3g %10.3g " % (Ncov[1,0],Ncov[2,1],Ncov[2,0]))
            except np.linalg.linalg.LinAlgError:
                mcpars = (mcstats['temperature']['mean'],mcstats['temperature']['standard deviation'],
                    mcstats['beta']['mean'],mcstats['beta']['standard deviation'],
                    mcstats['scale']['mean'],mcstats['scale']['standard deviation'])
                print("MCMC  %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g " % tuple(mcpars))
            except ValueError:
                mcpars = (mcstats['temperature']['mean'],mcstats['temperature']['standard deviation'],
                    mcstats['beta']['mean'],mcstats['beta']['standard deviation'],
                    mcstats['scale']['mean'],mcstats['scale']['standard deviation'])
                print("MCMC  %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g " % tuple(mcpars))


            pylab.figure(ii)
            pylab.clf()
            agpy.pymc_plotting.hist2d(bbmc,'temperature','beta',fignum=ii,varslice=(0,None,None),bins=20)
            ax=pylab.gca()
            ax.add_artist(matplotlib.patches.Ellipse(mp.params[:2], mp.perror[0], mp.perror[1],
                    mp.covar[1,0]/mp.covar[0,0]*90, edgecolor='green',
                    facecolor='none'))
            pylab.plot(temperature,beta,'kx')
            
            pylab.figure(ii+10)
            pylab.clf()

            input_flux = modified_blackbody_wavelength(wavelength, temperature,
                    beta=beta, wavelength_units='microns', normalize=False, logN=column,
                    logscale=16)
            recovered_flux = modified_blackbody_wavelength(wavelength, mcpars[0],
                    beta=mcpars[2], wavelength_units='microns', normalize=False, logN=mcpars[4],
                    logscale=16)

            pylab.plot(wavelength, flux, marker='o', linestyle='none')
            pylab.errorbar(wavelength, input_flux, err, linestyle='none')
            pylab.plot(wavelength, recovered_flux)

    for kk,errlevel in enumerate(errlevels):
        pylab.figure(100+kk)
        pylab.clf()
        for ii,(wavelength) in enumerate(wavelengths):
            flux = modified_blackbody_wavelength(wavelength, temperature,
                    beta=beta, wavelength_units='microns', normalize=False, logN=column,
                    logscale=16)
            err = flux*errlevel
            err[wavelength<1100] = flux[wavelength<1100]*0.02

            bbmc = fit_blackbody_montecarlo(wavelength, flux,
                    blackbody_function=modified_blackbody_wavelength,
                    return_MC=True, wavelength_units='microns', nsamples=1,
                    scale_guess=nguess, beta_guess=bguess, temperature_guess=tguess,
                    scale_keyword='logN', max_scale=30, err=err, logscale=16,
                    burn=0)
            mps = []
            pylab.figure(kk)
            pylab.clf()
            pylab.plot(temperature,beta,'kx')
            betas,betaerr,temps,temperr = [],[],[],[]
            for jj in xrange(5):
                flux = bbmc.flux.rand()
                mp = fit_blackbody(wavelength, flux, err=err,
                        blackbody_function=modified_blackbody_wavelength,
                        logscale=16, guesses=(tguess, bguess, nguess), wavelength_units='microns')
                
                mps.append(mp)
                
                temps.append(mp.params[0])
                betas.append(mp.params[1])
                betaerr.append(mp.perror[1])
                temperr.append(mp.perror[0])
            pylab.errorbar(temps,betas,xerr=temperr, yerr=betaerr, linestyle='none')
            pylab.figure(100+kk)
            print "%s sn=%f  beta=%f+/-%f" % (wavelength[-1],1/errlevel,np.mean(betas),np.std(betas))
            pylab.hist(betas,alpha=0.5,label="Longest Wavelength %s $\\mu m$ S/N=%0.1f" % (wavelength[-1],1/errlevel),histtype='stepfilled',bins=20)
            if ii==0:
                pylab.vlines(beta,*pylab.gca().get_ylim(),linestyle='--',color='k',label="Input value $\\beta=%f$" % beta)
            pylab.legend(loc='best')
            #pylab.savefig("/Users/adam/agpy/tests/longwav%i_sn%i_Herschelsn50_bb_test.png" % (wavelength[-1],1/errlevel))
