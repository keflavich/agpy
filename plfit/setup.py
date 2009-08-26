from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    ext_modules = [ 
      Extension("cplfit", ["cplfit.pyx"],
        include_dirs = [numpy.get_include(),
            '/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/Cython/',
            '.'],
        extra_compile_args=['-O3']
        )
      ],
    cmdclass = {'build_ext': build_ext},
)

