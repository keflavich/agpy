"""
===============
Dust emissivity
===============

nu is in GHz everywhere

"""
import blackbody,constants
from math import exp,log

def kappa(nu, nu0=constants.c/500e-4, kappa0=0.005, beta=1.75):
    return kappa0*(nu*1e9/nu0)**(-beta)

def snu(nu, column, kappa, temperature):
    tau = kappa / constants.mh
    snu = blackbody.blackbody(nu,temperature, normalize=False)
    return snu

def snudnu(nu, column, kappa, temperature, bandwidth):
    return snu(nu, column, kappa, temperature) * bandwidth

#def mass(nu, column, kappa, temperature, beamomega, distance=1):
#    return snu(nu, column, kappa, temperature) * beamomega / constants.msun * 1e23 / constants.mh * (distance*1000*constants.pc)**2

def snuofmass(nu, mass, beamomega, distance=1, temperature=20):
    column = mass * constants.msun / (beamomega * (distance*constants.kpc)**2)# g cm^-2 
    tau = kappa(nu) * column * beamomega
    bnu = blackbody.blackbody(nu, temperature, normalize=False, frequency_units='GHz')
    snu = bnu * (1.0-exp(-tau)) * 1e23
    return snu

def colofsnu(nu, snu, beamomega, temperature=20):
    bnu = blackbody.blackbody(nu, temperature, normalize=False, frequency_units='GHz')
    tau = -log(1-snu*1e-23 / bnu)
    column = tau / kappa(nu) / constants.mh / constants.muh2 / beamomega
    return column

def massofsnu(nu, snu, beamomega, distance=1, temperature=20):
    col = colofsnu(nu, snu, beamomega, temperature)
    mass = col * constants.mh * constants.muh2 * beamomega * (distance*constants.kpc)**2 
    return mass / constants.msun
