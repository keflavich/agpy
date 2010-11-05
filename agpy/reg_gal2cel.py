try:
    import coords
    import pyregion
except ImportError:
    print "reg_gal2cel requires coords & pyregion"
import posang

def gal2cel(regfile):
    """
    Converts a region file from galactic to celestial coordinates including
    position angle reference from the center of the box (right now only works
    on box regions)

    Requires pyregion with the ShapeList.write() function implemented...
    not clear if that exists in 1.0
    """
    reg = pyregion.open(regfile)

    for R in reg:
        if R.name == 'box':
            x,y,dx,dy,angle = R.coord_list

            posn = coords.Position([x,y],system='galactic')
            ra,dec = posn.j2000()

            newang = posang.posang(x-dx,y,x+dx,y,system='galactic')

            coord_list = [ra,dec,dx,dy,angle-newang-90]

            R.coord_format = 'fk5'
            R.coord_list = coord_list
            R.params = coord_list

    reg.write(regfile[:-4] + "_fk5.reg")

