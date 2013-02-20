import agpy
from agpy.blackbody import *
import numpy as np


wavelength = np.array([20,70,160,250,350,500,850,1100])
flux = modified_blackbody_wavelength(wavelength, 15, beta=1.75,
                wavelength_units='microns', normalize=False, logN=22, logscale=16)
err = 0.1 * flux
flux += np.random.randn(len(wavelength)) * err
tguess, bguess, nguess = 20.,2.,21.5
lm = fit_blackbody_lmfit(wavelength, flux, err=err,
                 blackbody_function=modified_blackbody_wavelength, logscale=16,
                 guesses=(tguess,bguess,nguess),
                 wavelength_units='microns')
print lm.params
       
# If you want to fit for a fixed beta, do this:
#import lmfit
parameters = lmfit.Parameters(OrderedDict([ (n,lmfit.Parameter(name=n,value=x)) for n,x
                in zip(('T','beta','N'),(20.,2.,21.5)) ]))

parameters['beta'].vary = False
lm = fit_blackbody_lmfit(wavelength, flux, err=err,
                 blackbody_function=modified_blackbody_wavelength, logscale=16,
                 guesses=parameters,
                 wavelength_units='microns')
print lm.params
