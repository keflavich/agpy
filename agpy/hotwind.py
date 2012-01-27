# formulae for hot wind computations
#
# Panagia & Felli 1975 / Wright & Barlow 1975:
from agpy.constants import *
from numpy import log10,exp,log

def mdot(snu, nu=10, Te=1e4, vinf=2000, mue=1.3, d=1):
    """
    [snu] = mJy
    nu = 10 GHz (frequency)
    Te = 10^4 K (ionized gas temp)
    vinf = 2000 km/s (maximum wind speed)
    mue = 1.3 (mean atomic weight per free electron)
    d = distance (kpc)

    return: mdot in Msun/yr
    """
    mdot = ( snu / 7.26 * (nu/10.)**-0.6 * (Te/1e4)**-0.1 * d**2 * (mue*vinf/100)**(4/3.) ) **(3/4.) * 1e-6
    return mdot

def mdotvinfr(snu, nu=10, Te=1e4, vinf=2000, mue=1.3, d=1, R=25):
    """
    [snu] = mJy
    nu = 10 GHz (frequency)
    Te = 10^4 K (ionized gas temp)
    vinf = 2000 km/s (maximum wind speed)
    mue = 1.3 (mean atomic weight per free electron)
    d = distance (kpc)

    return: Dmom in g*cm/s^2

    Sp. type Log D0 x alphaP
    A I 14.22 	 2.41 2.64 	 0.47 0.38 	 0.07
    Mid B I 17.07 	 1.05 1.95 	 0.20 0.51 	 0.05
    Early B I 21.24 	 1.38 1.34 	 0.25 0.75 	 0.15
    O I 20.69 	 1.04 1.51 	 0.18 0.66 	 0.06
    O III, V 19.87 	 1.21 1.57 	 0.21 0.64 	 0.06

    """
    mdot = ( snu / 7.26 * (nu/10.)**-0.6 * (Te/1e4)**-0.1 * d**2 * (mue*vinf/100)**(4/3.) ) **(3/4.) * 1e-6
    return mdot*msun/yr * (vinf*1e5) * (R)**0.5

def LofMdot(logdmom, x=1.51, D0=20.69, alphaP=0.66):
    """
    Kudritzki 2000
    log Dmom - log D0 = x log(L/Lsun); Dmom = Mdot * vinf *R/rsun

    returns L in log10(Lsun)
    """

    L = ( (logdmom-D0) / x )
    return L

def SnuOfMdot(nu, dkpc, mdot=1e-5, mu=1.2, vwind=1e3, Zbar=1, Te=1e4): 
    """
    Panagia and Felli 1975 equation 24

    [nu] = GHz
    [Te] = K
    [mdot] = msun/yr
    [vwind] = km/s
    mu = mean particle mass in AMU (presumably)
    Zbar = average ionic charge (assumed 1, but not sure that's good)
    dkpc = distance in kpc
    """

    snu = 5.12 * ( (nu/10.)**0.6 * (Te/1e4)**0.1 * (mdot/1e-5)**(4/3.) *
            (mu/1.2)**(-4/3.) * (vwind/1e3)**(-4/3.) * (Zbar)**(2/3.) *
            (dkpc)**-2.)

    # 1 mfu = 1 milli-flux-unit = 10^-29 W m^-2 Hz ^-1 = 1 mJy
    return snu
