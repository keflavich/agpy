#from distutils.core import setup
#from distutils.extension import Extension
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
#from numpy.distutils.core import build_ext
from Cython.Distutils import build_ext
import numpy

print "To create cplfit.so (for importing), call command: "
print "python setup.py build_ext --inplace"

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
    numpy_include_dirs = get_numpy_include_dirs()
except AttributeError:
    numpy_include_dirs = numpy.get_include()


ext_cplfit = Extension("cplfit", ["cplfit.pyx"], include_dirs = [numpy_include_dirs,
   '/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/Cython/','.'],
   extra_compile_args=['-O3'])
#ext_fplfit = Extension(name="fplfit",
#                    sources=["fplfit.f"])

if __name__=="main":
    setup(
        ext_modules = [ ext_cplfit ],
        cmdclass = {'build_ext': build_ext}
    )

print "I can't get numpy.distutils to compile the fortran.  To do it yourself, run some variant of:"
print 'f2py -c fplfit.f -m fplfit'

# try:
#     os.system('f2py -c fplfit.f -m fplfit')
# except:
#     print "Could not build fplfit"

