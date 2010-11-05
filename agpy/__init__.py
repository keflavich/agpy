"""
This is an init file.

The functions included below are the 'mature' codes from the agpy package.
"""
#from luminosity import luminosity
import readcol as readcol_mod
from readcol import readcol
#from UCHIIfitter import HIIregion
import gaussfitter
from gaussfitter import moments,twodgaussian,gaussfit,onedgaussian,onedgaussfit
import kdist
from kdist import kdist,vector_kdist
#from plfit import plfit
from reg_gal2cel import gal2cel
from posang import posang
import densitymap
import downsample as downsample_mod
from downsample import downsample,downsample_cube
import correlate2d as correlate2d_mod
from correlate2d import correlate2d
from psds import PSD2
import convolve as convolve_mod
from convolve import convolve,smooth
#from asinh_norm import AsinhNorm
from radialprofile import azimuthalAverage
import showspec
from contributed import parallel_map


# import all of the functions but not the modules...
__all__ = ['readcol','gaussfitter','kdist','reg_gal2cel','posang','densitymap','downsample','correlate2d','psds','convolve','radialprofile',
        'gal2cel',
        'convolve','smooth',
        'azimuthalAverage',
        'kdist','vector_kdist',
        'moments','twodgaussian','gaussfit','onedgaussian','onedgaussfit']
