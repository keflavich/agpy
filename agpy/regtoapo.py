try:
    import pyregion
    import coords
except ImportError:
    print "regtoapo requires pyregion and coords packages"

from region_positions import *

def regtoapo(regfile,outfile,filtercolor=None):
    """
    Take an input ds9 .reg file and create an output file
    formatted to work with APO's TUI

    filtercolor - if specified, will ignore any regions of this color
    """

    reg = pyregion.open(regfile)
    outf = open(outfile,'w')

    for r in reg:
        if r.attr[1]['color'] == filtercolor or not r.attr[1].has_key('text'):
            continue
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



if __name__ == '__main__':
    """
    Allows script to be run from the command line:
    python regtoapo.py blah.reg bleh.txt purple
    """
    import sys
    regtoapo(*sys.argv[1:])
   
