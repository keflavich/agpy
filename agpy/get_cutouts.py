"""
Make cutouts with all available data of a given position
"""
import glob
import pyfits
import pywcs
import coords
import os
import agpy.cubes,agpy.grep,agpy.cutout
import timer

@timer.print_timing
def findimages(xc, yc, searchpath, coordsys='celestial', verbose=False,
        ignore_imagetypes=('_area','power','angle','weight','residual','smooth','model','mask','noise','label','nhits','covg','std','rms')):
    """
    Given a single coordinate pair, searches for images that contain the coordinate
    """
    if verbose: print "Searching %s" % searchpath
    imlist = glob.glob(searchpath+"/*.fits")
    for ign in ignore_imagetypes: # remove undesired image types
        imlist = agpy.grep.grepv(ign, imlist)
    matchlist = []
    for fn in imlist:
        try: 
            head = pyfits.getheader(fn)
        except: # generic exception; there are lots of ways file reading could fail but I don't care
            continue
        head = agpy.cubes.flatten_header( head )
        if coords_in_image(xc,yc,head,coordsys=coordsys):
            matchlist.append(fn)
            if verbose:
                print "Matched %s" % fn

    if len(matchlist) > 1 and type(matchlist[0]) is list:
        # remove blanks / make into a single list
        return reduce(lambda a,b: a+b,matchlist,[])
    else:
        return matchlist

#@timer.print_timing # this gets called too many times
def coords_in_image(xc,yc,header,coordsys='celestial'):
    """
    Determine whether the coordinates are within the boundaries of the image
    """
    try:
        wcs = pywcs.WCS(header)
    except: # if the header isn't WCS compatible, we don't want it
        return False

    try:
        if coordsys=='celestial' and wcs.wcs.lngtyp=='GLON':
            xc,yc = coords.Position((xc,yc),system=coordsys).galactic()
        elif coordsys=='galactic' and wcs.wcs.lngtyp=='RA':
            xc,yc = coords.Position((xc,yc),system=coordsys).j2000()

        xp,yp = wcs.wcs_sky2pix(xc,yc,0)
    except:
        return False

    if xp > wcs.naxis1 or xp < 0:
        return False
    elif yp > wcs.naxis2 or yp < 0:
        return False
    else:
        return True

@timer.print_timing
def find_all_images(xc,yc,dirlist, flatten=False, **kwargs):
    """
    Given a list of directories, search all the directories for overlap with the coordinate
    """
    imlist = [findimages(xc,yc,dire,**kwargs) for dire in dirlist if dire is not []]
    if flatten and type(imlist[0]) is list:
        return reduce(lambda a,b: a+b, imlist, [])
    else:
        return imlist

@timer.print_timing
def get_cutouts(xcoord,ycoord,xwidth,ywidth, coordsys='galactic',
        ignore_imagetypes=('_area','power','angle','weight','residual','smooth','model','mask','noise','label','nhits','covg','std','rms'),
        flist='find', savedir=None, clobber=True, verbose=False, 
        **kwargs):
    """
    Create cutouts from all possible images in the searched directories.
    """

    if flist == 'find':
        # Get list
        flist = find_all_images(xcoord, ycoord, dirs, coordsys=coordsys,
                ignore_imagetypes=ignore_imagetypes, verbose=verbose,
                flatten=True, **kwargs)
    elif type(flist) not in (tuple,list):
        raise TypeError("flist must be a list (or tuple) of filenames")

    # filter out image types I don't want
    for ign in ignore_imagetypes:
        flist = agpy.grep.grepv(ign, flist)

    # make all the cutouts
    cutouts = []
    for fn in flist:
        try:
            co = agpy.cutout.cutout(fn, xcoord, ycoord, xwidth, ywidth, units='wcs', coordsys=coordsys, verbose=verbose)
        except agpy.cutout.DimensionError:
            header = pyfits.getheader(fn)
            wcs = pywcs.WCS( agpy.cubes.flatten_header(header) )
            xc,yc,xw,yw = agpy.cubes.aper_world2pix((xcoord,ycoord,xwidth,ywidth),wcs,coordsys=coordsys,wunit='degree')
            try:
                co = agpy.cubes.subcube( pyfits.getdata(fn), xcoord, xwidth, ycoord, ywidth, header=header, return_HDU=True, widthunits='pixels' )
            except Exception as ex:
                print "%s failed: " % fn, ex
                continue
        except Exception as ex:
            print "%s failed: " % fn, ex
            flist.remove(fn)
            continue
        cutouts.append(co)
        if savedir is not None:
            try:
                co.writeto("%s/%s" % (savedir, os.path.split(fn)[1].replace(".fits","_cutout.fits")),clobber=clobber)
            except pyfits.VerifyError:
                print "Could not verify FITS header for %s.  Image dimensions were %s" % (fn, co.data.shape)

    return cutouts

@timer.print_timing
def find_directories(rootdir, ignoredirs=('v0.7','powerspectra','AGBT','Frame')):
    """
    Find directories containing FITS files
    """

    paths = []
    for root,dirs,files in os.walk(str(rootdir)):
        for igd in ignoredirs:
            if igd in root:
                continue
        if ".fits" in (str(os.path.splitext(f)[1]) for f in files): 
            paths.append(root)

    return paths

default_rootdirs = ['/Volumes/disk2/data/','/Volumes/WD_2/HiGal/','/Volumes/disk4/']

#dirs = reduce(lambda a,b: a+b, [find_directories(dr) for dr in default_rootdirs], []) 

standard_dirs = [
        '/Volumes/disk2/data/MGPS/',
        '/Volumes/disk2/data/MSX/',
        '/Volumes/disk2/data/WISE/',
        '/Volumes/disk2/data/bally_CO/',
        '/Volumes/disk2/data/bgps/releases/IPAC/',
        '/Volumes/disk2/data/bgps/releases/v2.0/',
        '/Volumes/disk2/data/bgps/releases/v2.0/August2011/',
        '/Volumes/disk2/data/c2d/',
        '/Volumes/disk2/data/cara/csfiles/',
        '/Volumes/disk2/data/cara/glimpsev3/',
        '/Volumes/disk2/data/co/',
        '/Volumes/disk2/data/glimpse/',
        '/Volumes/disk2/data/glimpseii/',
        '/Volumes/disk2/data/grs/',
        '/Volumes/disk2/data/harp/fits/',
        '/Volumes/disk2/data/iras/',
        '/Volumes/disk2/data/magpis/',
        '/Volumes/disk2/data/mips/',
        '/Volumes/disk2/data/mipsgal/',
        '/Volumes/disk2/data/mosaics/',
        '/Volumes/disk2/data/motte/',
        '/Volumes/disk2/data/nvss/',
        '/Volumes/disk2/data/scuba/',
        '/Volumes/disk2/data/sharc/processd/',
        '/Volumes/disk2/data/vgps/',
        '/Volumes/WD_2/HiGal/',
        '/Volumes/disk4/cmz/',
        '/Volumes/disk4/gc/',
        '/Volumes/disk4/higal-gc/',
        '/Volumes/disk4/l44_kirk/',
        '/Volumes/disk4/mosaics/',
        '/Volumes/disk4/nvss/',
        '/Volumes/disk4/orion/',
        '/Volumes/disk4/perseus/',
        '/Volumes/disk4/vlss/',
        ]

dirs=standard_dirs
