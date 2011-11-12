try:
    import pyregion
    import coords
except ImportError:
    print "regtocima requires pyregion and coords packages"
import numpy

from region_positions import *

def regtocima(regfile,outfile,filtercolor=None):
    """
    Take an input ds9 .reg file and create an output file
    formatted to work with Arecibo's CIMA
    [NOT FUNCTIONAL]

    filtercolor - if specified, will ignore any regions of this color
    """

    print "Doesn't work yet"
    return

    reg = pyregion.open(regfile)
    outf = open(outfile,'w')

    for r in reg:
        if r.attr[1]['color'] == filtercolor or not r.attr[1].has_key('text'):
            continue

        text = pos_to_name(r)

        if r.name == 'box':
            x,y,dx,dy,posang = r.coord_list
            radec = position_region(r).hmsdms()
            # ds9's position angles are 90 degrees offset from APO's
            if posang+90 > 360: posang -= 360
            print >>outf,'%32s %26s rotangle=%f' % ( '"%s"' % r.attr[1]['text'],radec,posang+90)
        else:
            radec = position_region(r).hmsdms()
            print >>outf,'%32s %26s' % ('"%s"' % r.attr[1]['text'],radec)

    outf.close()

def regtouptime(regfile,outfile,maxza=15.0):
    """
    Converts a region file to a list of target names + how long they are observable @ arecibo

    maxza = 15.0 degrees for "broken" Arecibo (2010) or 18-19 for normal Arecibo
    """

    # coordinates @ http://www.naic.edu/~astro/guide/node2.html
    ra_arecibo  = (66+45/60.+11.1/3600.)
    dec_arecibo = (18+20/60.+36.6/3600.)

    reg = pyregion.open(regfile)
    outf = open(outfile,'w')

    for r in reg:
        #if r.attr[1]['color'] == filtercolor or not r.attr[1].has_key('text'):
        #    continue
        text = pos_to_name(r)
        radec = position_region(r).hmsdms()
        ra,dec = position_region(r).j2000()
        za_min = numpy.abs(dec-dec_arecibo)
        # isoscoles triangle with za_min as perpendicular bisector, maxza as long leg lenth
        tup = numpy.sqrt(maxza**2-za_min**2)*2.0/15.0 * 60.0 # minutes
        print >>outf,'%12s %10.2f %10.1f' % (text,za_min,tup)

    outf.close()

if __name__ == '__main__':
    """
    Allows script to be run from the command line:
    python regtocima.py blah.reg bleh.txt purple
    """
    import sys
    if len(sys.argv) > 1:
        regtocima(*sys.argv[1:])
   
