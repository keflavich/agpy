import numpy
try:
    import coords

    def position_region(reg):
        """
        small wrapper to get a Position object using the correct coordinate system

        reg must by a pyregion Shape instance
        """
        x,y = reg.coord_list[:2]
        posn = coords.Position([x,y],system=coords_format(reg.coord_format))
        return posn

    def coords_format(format):
        """
        Convert from ds9's 'fk5' and 'icrs' naming convention to the
        'celestial'/'galactic' used by coords
        """
        if format == 'galactic':
            return 'galactic'
        elif format in ['fk5','icrs']:
            return 'celestial'

    def pos_to_name(reg):
        """
        Given a region, returns a name based on Galactic coordinates
        """
        l,b = position_region(reg).galactic()           
        if numpy.sign(b) == 1:
            pm = "+"
        else:
            pm = "-"
        text = "G%4.2f%1s%4.2f" % (l,pm,abs(b))
        return text

except ImportError:
    print "region_positions requires coords"

