#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
import pyregion
import pyfits
import ds9
import numpy

def ds9_photometry(xpapoint):
    D = ds9.ds9(xpapoint)
    reg = pyregion.parse(D.get("regions selected -format ds9 -system wcs -sky fk5 -skyformat sexagesimal"))
    pf = D.get_pyfits()
    mask = reg.get_mask(pf[0])
    arr = pf[0].data
    return arr[mask].sum(),arr[mask].mean(),numpy.median(arr[mask]),arr[mask].std(),mask.sum()

if __name__ == "__main__":

    import sys
    xpaname = sys.argv[1]

    print "Sum: %g  Mean: %g  Median: %g  RMS: %g  NPIX: %i" % ds9_photometry(xpaname)

