try:
    import pyregion
    import coords
except ImportError:
    print "regtoapo requires pyregion and coords"

def regtoapo(regfile,outfile,filtercolor=None):

    reg = pyregion.open(regfile)
    outf = open(outfile,'w')

    for r in reg:
        if r.attr[1]['color'] == filtercolor or not r.attr[1].has_key('text'):
            continue
        if r.name == 'box':
            x,y,dx,dy,posang = r.coord_list
            radec = position_region(r).hmsdms()
            print >>outf,'"%32s" %26s rotangle=%f' % (r.attr[1]['text'],radec,posang)
        else:
            radec = position_region(r).hmsdms()
            print >>outf,'"%32s" %26s' % (r.attr[1]['text'],radec)

    outf.close()


def position_region(reg):                
    x,y = reg.coord_list[:2]
    posn = coords.Position([x,y],system=coords_format(reg.coord_format))
    return posn

def coords_format(format):
    if format == 'galactic':
        return 'galactic'
    elif format == 'fk5':
        return 'celestial'

if __name__ == '__main__':
    import sys
    regtoapo(*sys.argv[1:])
   
