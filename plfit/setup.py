#from distutils.core import setup
#from distutils.extension import Extension
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from Cython.Distutils import build_ext
import numpy

print "To create cplfit.so (for importing), call command: "
print "python setup.py build_ext --inplace"

ext_cplfit = Extension("cplfit", ["cplfit.pyx"], include_dirs = [numpy.get_include(),
   '/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/Cython/','.'],
   extra_compile_args=['-O3'])
ext_fplfit = Extension(name="fplfit",
                    sources=["fplfit.f"])

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

