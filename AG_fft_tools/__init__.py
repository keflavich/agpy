"""
=========
FFT Tools
=========

Wrappers around numpy, scipy, and pyfftw tools to perform 2D convolution in
general, smoothing with a set of 'standard' kernels, and computing power
spectra and PSDs.

"""
import warnings
warnings.warn("agpy.AG_fft_tools is deprecated."
              " Instead, seek out image_tools (which contains all these tools): "
              "https://github.com/keflavich/image_tools",
              DeprecationWarning)

from correlate2d import correlate2d
from psds import PSD2
from smooth_tools import smooth
from convolve_nd import convolvend
from convolve_nd import convolvend as convolve
import fast_ffts
from upsample import dftups,upsample_image
from shift import shift
