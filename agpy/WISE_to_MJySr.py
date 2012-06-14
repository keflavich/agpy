"""
Small script to convert a WISE image from DN to MJy/sr
"""
import pyfits
import numpy as np

dn_to_jy = {
    1:1.9350E-06,
    2:2.7048E-06,
    3:1.8326e-06,
    4:5.2269E-05,
    }

beam_area = {
    1:  6.08   * 5.60 * np.pi / (8*np.log(2)) / 206265.**2,
    2:  6.84   * 6.12 * np.pi / (8*np.log(2)) / 206265.**2,
    3:  7.36   * 6.08 * np.pi / (8*np.log(2)) / 206265.**2,
    4:  11.99  * 11.65* np.pi / (8*np.log(2)) / 206265.**2,
    }
        
bmaj_min_pa = {
    1:(6.08,5.60,3),
    2:(6.84,6.12,15),
    3:(7.36,6.08,6),
    4:(11.99,11.65,0),
    }

        

def WISE_to_MJySr(filename, outfilename=None, clobber=False):
    """
    Convert a WISE image from DN to MJy/Sr as per
    http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_3f.html

    *WARNING* The WISE documentation says that the images cannot be calibrated
    to surface brightness!  Flux measurements extracted from the WISE data are
    likely to be incorrect as a result of these modifications!  But, for
    qualitative purposes, it should be fine.

    Parameters
    ----------
    filename - string
        The name of the WISE .fits file acquired from IPAC
        (must have keywords BUNIT='DN' and MAGZP)
    outfilename - string (optional)
        If specified, will write to this file name.  Otherwise,
        will try to overwrite input filenmae
    clobber - bool
        Overwrite specified output file name?
    """
    fitsfile = pyfits.open(filename)
    header = fitsfile[0].header
    band = header['BAND']

    if header['BUNIT'].strip() != 'DN':
        raise TypeError("Wrong BUNIT %s" % header['BUNIT'])

    if band not in dn_to_jy:
        raise ValueError("WISE only has 4 bands!  This file has band=%s" % band)

    fitsfile[0].data *= dn_to_jy[band] * 1e-6 / beam_area[band]
    header['BUNIT'] = 'MJy/sr'
    header.update('OMEGABM',beam_area[band])
    header.update('DNTOJY',dn_to_jy[band])
    bmaj,bmin,bpa = bmaj_min_pa[band]
    header.update('BMAJ',bmaj)
    header.update('BMIN',bmin)
    header.update('BPA', bpa)

    if outfilename is not None:
        fitsfile.writeto(outfilename, clobber=clobber)
    else:
        fitsfile.writeto(filename, clobber=clobber)

if __name__ == "__main__":
    import optparse

    parser=optparse.OptionParser()

    parser.add_option("--clobber",help="Overwrite output file if it exists?",action='store_true',default=False)
    parser.add_option("--outfilename","-o",help="Output file name (optional)",default=None)
    parser.set_usage("%prog WISE_file.fits [options]")
    parser.set_description(
    """
    Convert a WISE image in DN to MJy/Sr 
    """)

    options,args = parser.parse_args()

    if len(args) > 0 and options.outfilename is None and options.clobber:
        for fn in args:
            WISE_to_MJySr(fn, clobber=options.clobber)
    elif len(args) > 0:
        raise ValueError("Too many filenames; can only use multiple file names if clobber=True and no outfilename specified.")
    else:
        WISE_to_MJySr(*args, clobber=options.clobber)
