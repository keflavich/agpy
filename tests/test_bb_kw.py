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


