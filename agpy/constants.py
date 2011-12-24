"""
=========
Constants
=========

A simple list of important constants in CGS units

"""

#constants
hplanck = h = 6.626068e-27    # Planck's Constant (erg s)
k = kb = kB = 1.3806503e-16   # Boltzmann's Constant (erg / K)
G = 6.67300e-8                # Newton's constant (cm^3 g^-1 s^-2)
c = 2.99792458e10             # Speed of Light (cm/s)
m_co = 4.65e-23               # Mass of CO (g)
m_h2 = 3.32e-24               # Mass of H_2 (g)
msun = 1.9889e33              # grams
lsun = 3.839e33               # erg/s
h2toCO = 1e4                  # X_CO = CO abundance
kmstocms  = 1e5               # conversion from km/s to cm/s
au = 1.496e13                 # 1 astronomical unit in cm
pc = 3.086e18                 # 1 parsec in cm
kpc = pc * 1e3
muh=1.4                       # Mass conversion factor for H 
muh2=2.8                      # Mass conversion factor for H2 (see Kauffman 2009)
yr = year = 86400*365.        # 1 year in seconds
amu = mh = 1.66053886e-24     # 1 atomic mass unit (1 hydrogen atom) in grams
alphab = 2.6e-13              # Case B recombination coefficient for H
pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679 # why not?
#                                            ^
# According to wikipedia, we only need accuracy to the 39th digit to make a
# circle the size of the observable universe accurate to the size of a hydrogen
# atom                                           
# Hmm... r_H ~ 0.5 angstroms ~ 5*10^-9 cm
#        r_Uni ~ 10 Gpc (just to be round...) ~ 10^10 * 3*10^18 = 3*10^28 cm
#        r_Uni/r_H ~ 10^37
#        Close enough.

import numpy as np
FWHM = np.sqrt(8*np.log(2))
HWHM = FWHM/2.
