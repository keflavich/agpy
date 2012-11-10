#!/usr/bin/env python -W ignore::DeprecationWarning
#!/Library/Frameworks/Python.framework/Versions/2.6/bin/python
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pyfits
    import pywcs
import pyregion
import ds9
import numpy
import sys
import re

def ds9_photometry(xpapoint):
    D = ds9.ds9(xpapoint)
    try:
        reg = pyregion.parse(D.get("regions selected -format ds9 -system wcs -sky fk5 -skyformat sexagesimal"))
    except Exception as ex:
        print ex
        raise ex
    pf = D.get_pyfits()
    ph = pyfits.PrimaryHDU(data=pf[0].data,header=pf[0].header)
    mask = reg.get_mask(ph)
    arr = pf[0].data
    wherenotnan = (arr == arr)
    mask = mask*wherenotnan
    hdr = pf[0].header
    try:
        wcs = pywcs.WCS(hdr.tostring())
    except AttributeError:
        wcs = pywcs.WCS(hdr)
    try:
        try:
            bmaj = float(hdr['BMAJ'])
            bmin = float(hdr['BMIN'])
        except KeyError:
            # VLA imfits
            bmin = None; bmaj = None
            for k,v in hdr.iteritems():
                if numpy.iterable(v) and "BMAJ" in v:
                    bmaj = float(v.split()[3])
                    bmin = float(v.split()[5])
            if bmin is None or bmaj is None:
                raise KeyError("BMIN and BMAJ not found")
        try:
            cd1 = wcs.wcs.cd[0,0]
            cd2 = wcs.wcs.cd[1,1]
        except AttributeError:
            cd1,cd2 = wcs.wcs.cdelt[:2]
        ppbeam = 2*numpy.pi*bmin*bmaj / abs(cd1*cd2) / (8*numpy.log(2))
        #print "CD1: %g  CD2: %g" % (cd1, cd2)
        sys.stdout.write( "BMAJ: %g  BMIN: %g  PPBEAM: %g   SUM/PPBEAM: %g\n" % (bmaj,bmin,ppbeam,arr[mask].sum()/ppbeam) )
    except KeyError:
        print "ds9_phot failed - check for BMAJ/BMIN in header"
        pass
    except Exception as inst:
        print "ds9_phot failed - not a header KeyError, something else"
        print inst.args
        pass
    return arr[mask].sum(),arr[mask].mean(),numpy.median(arr[mask]),arr[mask].std(),mask.sum()

if __name__ == "__main__":

    import sys
    xpaname = sys.argv[1]

    sys.stdout.write( "Sum: %g  Mean: %g  Median: %g  RMS: %g  NPIX: %i\n" % ds9_photometry(xpaname) )

    sys.exit()
