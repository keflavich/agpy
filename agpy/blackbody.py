"""
Simple black-body calculator
"""
try:
    from numpy import exp
except ImportError:
    from math import exp

unitdict = {
        'cgs':{ 'h':6.626068e-27,
            'k':1.3806503e-16,
            'c':2.99792458e10},
        'mks':{ 'h':6.626068e-34,
            'k':1.3806503e-23,
            'c':2.99792458e8}
        }

frequency_dict = {
        'Hz':1.0,
        'kHz':1e3,
        'MHz':1e6,
        'GHz':1e9,
        'THz':1e12,
        }

def blackbody(nu,temperature,units='cgs',frequency_units='Hz'):
    # load constants in desired units
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']

    # convert nu to Hz
    nu = nu * frequency_dict[frequency_units]

    I = 2*h*nu**3 / c**2 * (exp(h*nu/(k*temperature)) - 1)**-1

    return I

wavelength_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'angstroms':1e-10,'A':1e-10,'Angstroms':1e-10,
        }

def blackbody_wavelength(lam,temperature,units='cgs',wavelength_units='Angstroms'):
    # load constants in desired units
    h,k,c = unitdict[units]['h'],unitdict[units]['k'],unitdict[units]['c']

    # convert nu to Hz
    lam = lam * wavelength_dict[wavelength_units] / (1e-2 if units=='cgs' else 1)

    I = 2*h*c**2 / lam**5 * (exp(h*c/(k*temperature*lam)) - 1)**-1

    return I
