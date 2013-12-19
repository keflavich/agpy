
# stupid script to install python after breaking it (which is a dumb thing to do)
#

cd ~/repos/matplotlib.git
git pull
export CC=/Developer/usr/bin/gcc-4.2
export CXX=/Developer/usr/bin/g++-4.2
export CFLAGS=-arch x86_64
make -f make.osx PREFIX=/Users/adam/repos/mpl_dependencies/ PYVERSION=2.7 deps mpl_install

cd ~/repos
for setupdir in aipy aplpy.git asciitable astropysics-dev atpy coords-svn cython Django-1.3.1 guppy h5py h5py-1.3.1 healpy-0.9.10 hyperspy idlsave Imaging-1.1.7 ipython-0.10.2 kapteyn-2.1 line_profiler lineid_plot mercurial-2.0 mpi4py-read-only mpmath-read-only mutagen pandas PFits psi pupynere-1.0.15 pyavm pyds9-1.0 pyephem pyfft-0.3.6 PyFFTW3-0.2.1 pyfits-2.4.0 pyregion-github pyspeckit python-dateutil-1.5 python-montage.git python-quantities python-sao-read-only pywcs-1.10-4.7 pywt s3cmd-1.0.1 sphinx thg vo-0.6 
do
    cd $setupdir
    python setup.py install
    cd ..
done

