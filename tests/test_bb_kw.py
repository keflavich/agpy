import agpy
import numpy as np
wavelength = np.array([20,70,160,250,350,500,850,1100])
flux = agpy.blackbody.modified_blackbody_wavelength(wavelength, 15, beta=2,
                logN=22, wavelength_units='microns', normalize=False,
                logscale=16)

err = 0.1 * flux
flux += np.random.randn(len(wavelength)) * err
tguess, bguess, nguess = 20.,2.,21.5
mp = agpy.blackbody.fit_blackbody(wavelength, flux, err=err,
	blackbody_function=agpy.blackbody.modified_blackbody_wavelength,
	guesses=(tguess, bguess, nguess),logscale=16,
	wavelength_units='microns')
print mp.params
print mp.perror


# logscale=16, ???
# Temp, beta, n(H2)



wavelength = np.array([20,70,160,250,350,500,850,1100])
flux = agpy.blackbody.modified_blackbody_wavelength(wavelength, 15, beta=1.75,
                wavelength_units='microns', normalize=False, logN=22, logscale=16)
err = 0.1 * flux
flux += np.random.randn(len(wavelength)) * err
tguess, bguess, nguess = 20.,2.,21.5
lm = agpy.blackbody.fit_blackbody_lmfit(wavelength, flux, err=err,
                 blackbody_function=agpy.blackbody.modified_blackbody_wavelength, logscale=16,
                 guesses=(tguess,bguess,nguess),
                 wavelength_units='microns')
print lm.params
        

"""
Incomplete testing so far.  Looks like luminosity comes out very wrong.  Other
parameters seem to be recovered well enough.
"""
wavelength = np.array([20,70,160,250,350,500,850,1100])
frequency = 3e14 / wavelength
flux = agpy.blackbody.modified_blackbody(frequency, 15, beta=1.75,
                                         normalize=False, logN=22, logscale=16)
err = 0.1 * flux
flux += np.random.randn(len(wavelength)) * err
tguess, bguess, nguess = 20.,2.,21.5
mc = agpy.blackbody.fit_blackbody_montecarlo(frequency, flux, err=err,
                                             blackbody_function=agpy.blackbody.modified_blackbody,
                                             logscale=16,
                                             temperature_guess=tguess,
                                             beta_guess=bguess,
                                             scale_guess=nguess,
                                             #guesses=(tguess,bguess,nguess),
                                             scale_keyword='logscale',
                                             return_MC=True)
        
from agpy import pymc_plotting
pymc_plotting.hist2d(mc, mc.temperature, mc.beta)
