"""
====
agpy
====

The functions included below are the 'mature' codes from the agpy package.

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

"""
##from luminosity import luminosity
#import readcol as readcol_mod
from readcol import readcol
##from UCHIIfitter import HIIregion
import gaussfitter
from gaussfitter import moments,twodgaussian,gaussfit,onedgaussian,onedgaussfit
#import kdist
from kdist import kdist,vector_kdist
#from plfit import plfit
from reg_gal2cel import gal2cel
from posang import posang
#import densitymap
#import downsample as downsample_mod
from downsample import downsample,downsample_cube

#from asinh_norm import AsinhNorm
#import showspec # imports matplotlib = BAD
from contributed import parallel_map
from timer import print_timing
from region_photometry import region_photometry
from region_photometry_files import region_photometry_files
from PCA_tools import efuncs,pca_subtract,unpca_subtract,smooth_waterfall
import constants
import blackbody

import AG_fft_tools
from AG_fft_tools import *

import AG_image_tools
from AG_image_tools import *

import cutout
import get_cutouts

# import all of the functions but not the modules...
__all__ = ['readcol', 'gaussfitter', 'kdist', 'reg_gal2cel', 'posang',
    'densitymap', 'downsample', 'correlate2d', 'psds', 'convolve', 'radialprofile',
    'constants' 'gal2cel', 'convolve', 'smooth', 'azimuthalAverage',
    'azimuthalAverageBins', 'kdist', 'vector_kdist', 'moments', 'twodgaussian',
    'gaussfit', 'onedgaussian', 'onedgaussfit']
