#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
try:
    import pyregion
except ImportError:
    print "Region Photometry requires pyregion"
try:
    import pyfits
except ImportError:
    print "Region photometry requires pyfits"
import pyfits
import numpy
from agpy import mad

def region_photometry(regionfile,fitsfilename,outfile='/dev/tty',doprint=True):

    
    fitsfile = pyfits.open(fitsfilename)
    data = fitsfile[0].data
    gooddata = data==data
    header = fitsfile[0].header
    if header.has_key('BMAJ') and header.has_key('BMIN'):
        bmaj = float(header['BMAJ'])
        bmin = float(header['BMIN'])
        if header.has_key('CDELT1'):
            CD1 = header.get('CDELT1')
            CD2 = header.get('CDELT2')
        elif header.has_key('CD1_1'):
            CD1 = header.get('CD1_1')
            CD2 = header.get('CD2_2')
        else:
            CD1 = CD2 = 1
        ppbeam = 2*numpy.pi*bmin*bmaj / abs(CD1*CD2) / (8*numpy.log(2))

    reglist = pyregion.open(regionfile)

    if doprint:
        outf = open(outfile,'w')
        print >> outf," ".join(["%10s" % s for s in 
            ["Name","Total","Sum/ppbm","Mean","Median","StdDev","MAD-STD","NPIX"]])

    props = []

    for ii,reg in enumerate(reglist):
        regL = pyregion.ShapeList()
        regL.append(reg)
        regmask = regL.get_mask(hdu=fitsfile[0])*gooddata
        if regmask.sum() > 0:
            total,avg,med,std,npix = data[regmask].sum(),data[regmask].mean(),numpy.median(data[regmask]),data[regmask].std(),regmask.sum()
            rmad = mad.MAD(data[regmask])
        else:
            total,avg,med,std,rmad,npix = 0,0,0,0,0,0
        if reg.attr[1].has_key('text'):
            name = reg.attr[1]['text']
        else:
            name = "%i" % ii
        if doprint:
            print >> outf,"%10s" % (name), " ".join(["%10.5g" % s for s in 
                [total,total/ppbeam,avg,med,std,rmad,npix]])
        props.append([name,total,total/ppbeam,avg,med,std,rmad,npix])

    if doprint: outf.close()

    return props

if __name__ == "__main__":

    import sys
    regionfilename,fitsfilename = sys.argv[1:3]
    if len(sys.argv) > 4:
        outfilename = sys.argv[3]
    else:
        outfilename = "/dev/tty"

    region_photometry(regionfilename,fitsfilename,outfile=outfilename)

    sys.exit()
