"""
=========
FFT Tools
=========

Wrappers around numpy, scipy, and pyfftw tools to perform 2D convolution in
general, smoothing with a set of 'standard' kernels, and computing power
spectra and PSDs.

"""

from correlate2d import correlate2d
from psds import PSD2
from convolve import convolve,smooth
