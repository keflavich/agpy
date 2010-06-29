#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
import pyregion
import pyfits
import pywcs
import ds9
import numpy
import sys

def ds9_photometry(xpapoint):
    D = ds9.ds9(xpapoint)
    reg = pyregion.parse(D.get("regions selected -format ds9 -system wcs -sky fk5 -skyformat sexagesimal"))
    pf = D.get_pyfits()
    mask = reg.get_mask(pf[0])
    arr = pf[0].data
    wherenotnan = (arr == arr)
    mask = mask*wherenotnan
    hdr = pf[0].header
    wcs = pywcs.WCS(hdr)
    try:
        bmaj = float(hdr['BMAJ'])
        bmin = float(hdr['BMIN'])
        cdelt1,cdelt2 = wcs.wcs.cdelt[:2]
        cd1 = cdelt1 * wcs.wcs.cd[0,0]
        cd2 = cdelt2 * wcs.wcs.cd[1,1]
        ppbeam = 2*numpy.pi*bmin*bmaj / abs(cd1*cd2) / (8*numpy.log(2))
        # print "CD1: %g  CD2: %g" % (cd1, cd2)
        sys.stdout.write( "BMAJ: %g  BMIN: %g  PPBEAM: %g   SUM/PPBEAM: %g\n" % (bmaj,bmin,ppbeam,arr[mask].sum()/ppbeam) )
    except:
        print "ds9_phot failed - check for BMAJ/BMIN in header"
        pass
    return arr[mask].sum(),arr[mask].mean(),numpy.median(arr[mask]),arr[mask].std(),mask.sum()

if __name__ == "__main__":

    import sys
    xpaname = sys.argv[1]

    sys.stdout.write( "Sum: %g  Mean: %g  Median: %g  RMS: %g  NPIX: %i\n" % ds9_photometry(xpaname) )

