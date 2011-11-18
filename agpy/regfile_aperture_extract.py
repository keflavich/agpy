import pyregion
import pyfits
import numpy as np
import pywcs
from agpy.region_positions import *

def get_fluxes(regfile, outfile, inneraprad=35, outeraprad=60, hdu=None, PPBEAM=1.0, debug=False, print_nulls=False):
    """
    Extract fluxes from a region-defined aperture with inner and outer circular apertures
    specififed

    MUST BE IN GALACTIC COORDINATES
    """
    if hdu is None:
        raise ValueError("hdu keyword is required")

    data = hdu.data
    header = hdu.header
    wcs = pywcs.WCS(header)
    glonmax = wcs.wcs_pix2sky(0,0,0)[0]
    glonmin = wcs.wcs_pix2sky(data.shape[1],0,0)[1]

    reglist = pyregion.open(regfile)

    outf = open(outfile,'w')

    print "".join("%16s" % s for s in ['Source_Name','SumJy','ApSumJy','MeanApSumJy','SumJyBm','ApSumJyBm','BgMed','BgMean','BgStd','FracErrBg'])
    print >>outf,"".join("%16s" % s for s in ['Source_Name','SumJy','ApSumJy','MeanApSumJy','SumJyBm','ApSumJyBm','BgMed','BgMean','BgStd','FracErrBg'])

    for reg in reglist:
        glon,glat = position_region(reg).galactic()
        if not((glon > glonmin) and (glon < glonmax)):
            # these are the limits of the survey
            if print_nulls:
                print >>outf,"%16s" % sourcename,"".join("%16s" % f 
                        for f in ['-','-','-','-','-','-','-','-','-'])
            continue
        xc,yc = wcs.wcs_sky2pix(glon,glat,0)
        if xc < 0 or xc > data.shape[1] or yc < 0 or yc > data.shape[0]:
            if print_nulls:
                print >>outf,"%16s" % sourcename,"".join("%16s" % f 
                        for f in ['-','-','-','-','-','-','-','-','-'])
            continue
        regL = pyregion.ShapeList()
        reg.name = 'circle'
        while len(reg.coord_list) < 3:
            reg.coord_list.append(0)
        reg.coord_list[2] = inneraprad/3600.0  # set inner aperture (foreground) to R=25"
        regL.append(reg)
        innerap = regL.get_mask(hdu=hdu)
        if innerap.sum() == 0:
            print "Skipped a source that was in the boundaries: ",reg
            continue
        regL[0].coord_list[2] = outeraprad/3600.0  # set outer aperture (background) to R=100"
        #regL.append(reg) # I think this turns two circles into a panda?
        outerap = regL.get_mask(hdu=hdu)
        backreg = outerap-innerap

        total = data[innerap].sum()
        background = np.median(data[backreg])
        backmean = data[backreg].mean()
        backstd = data[backreg].std()

        sourcename = pos_to_name(reg)

        if backstd > total or total < 0:
            print "%s set to zero" % reg.attr[1]['text']
            total = 0
            total_backsub  = 0
            total_mbacksub = 0
        else:
            total_backsub  = total - innerap.sum() * background
            total_mbacksub = total - innerap.sum() * backmean

        print "%16s" % sourcename,"".join("%16.5g" % f 
                for f in [total/PPBEAM,total_backsub/PPBEAM,total_mbacksub/PPBEAM,total,total_backsub,background,backmean,backstd,backstd/total_backsub])
        print >>outf,"%16s" % sourcename,"".join("%16.5g" % f 
                for f in [total/PPBEAM,total_backsub/PPBEAM,total_mbacksub/PPBEAM,total,total_backsub,background,backmean,backstd,backstd/total_backsub])

    outf.close()

