"""
Small script to convert a MSX image from DN to MJy/sr
"""
import pyfits
import numpy as np

WMtoJy = {
    'A':7.133e6,
    'C':2.863e7,
    'D':3.216e7,
    'E':2.476e7,
}
        
wavelength_to_band = {
    8.28e-6 :'A',
    1.213e-5:'C',
    1.465e-5:'D',
    2.134e-5:'E',
    }

def MSX_to_MJySr(filename, outfilename=None, clobber=False, output_verify='fix'):
    """
    Convert an MSX image from DN to MJy/Sr as per
    http://irsa.ipac.caltech.edu/applications/MSX/MSX/imageDescriptions.htm

    Parameters
    ----------
    filename - string
        The name of the MSX .fits file acquired from IPAC
        (must have keywords BUNIT='W/m^2-sr' and WAVELENG)
    outfilename - string (optional)
        If specified, will write to this file name.  Otherwise,
        will try to overwrite input filenmae
    clobber - bool
        Overwrite specified output file name?
    output_verify - string
        pyfits' output_verify argument
    """
    fitsfile = pyfits.open(filename)
    header = fitsfile[0].header
    if header['BUNIT'].strip() != 'W/m^2-sr':
        raise TypeError("Wrong BUNIT %s" % header['BUNIT'])

    if header['WAVELENG'] not in wavelength_to_band:
        raise ValueError("MSX only has 4 bands!  This file has WAVELENG=%s" % header.get('WAVELENG'))

    band = wavelength_to_band[ header['WAVELENG'] ]

    fitsfile[0].data *= WMtoJy[band]
    header['BUNIT'] = 'MJy/sr'

    if outfilename is not None:
        fitsfile.writeto(outfilename, clobber=clobber, output_verify=output_verify)
    else:
        fitsfile.writeto(filename, clobber=clobber, output_verify=output_verify)

if __name__ == "__main__":
    import optparse

    parser=optparse.OptionParser()

    parser.add_option("--clobber",help="Overwrite output file if it exists?",action='store_true',default=False)
    parser.add_option("--outfilename","-o",help="Output file name (optional)",default=None)
    parser.set_usage("%prog MSX_file.fits [options]")
    parser.set_description(
    """
    Convert a MSX image in DN to MJy/Sr 
    """)

    options,args = parser.parse_args()

    if len(args) > 0 and options.outfilename is None and options.clobber:
        for fn in args:
            MSX_to_MJySr(fn, clobber=options.clobber)
    elif len(args) > 0:
        raise ValueError("Too many filenames; can only use multiple file names if clobber=True and no outfilename specified.")
    else:
        MSX_to_MJySr(*args, clobber=options.clobber)
