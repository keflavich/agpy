"""
Make cutouts with all available data of a given position
"""
import glob
import pyfits
import pywcs
import coords
import os
import agpy.cubes,agpy.grep,agpy.cutout

def findimages(xc,yc,searchpath,coordsys='celestial'):
    """
    Given a single coordinate pair, searches for images that contain the coordinate
    """
    imlist = glob.glob(searchpath+"/*.fits")
    matchlist = []
    for fn in imlist:
        try: 
            head = pyfits.getheader(fn)
        except: # generic exception; there are lots of ways file reading could fail but I don't care
            continue
        head = agpy.cubes.flatten_header( head )
        if coords_in_image(xc,yc,head,coordsys=coordsys):
            matchlist.append(fn)

    return matchlist

def coords_in_image(xc,yc,header,coordsys='celestial'):
    """
    Determine whether the coordinates are within the boundaries of the image
    """
    try:
        wcs = pywcs.WCS(header)
    except: # if the header isn't WCS compatible, we don't want it
        return False

    if coordsys=='celestial' and wcs.wcs.lngtyp=='GLON':
        xc,yc = coords.Position((xc,yc),system=coordsys).galactic()
    elif coordsys=='galactic' and wcs.wcs.lngtyp=='RA':
        xc,yc = coords.Position((xc,yc),system=coordsys).j2000()

    xp,yp = wcs.wcs_sky2pix(xc,yc,0)

    if xp > wcs.naxis1 or xp < 0:
        return False
    elif yp > wcs.naxis2 or yp < 0:
        return False
    else:
        return True

def find_all_images(xc,yc,dirlist, flatten=False, **kwargs):
    """
    Given a list of directories, search all the directories for overlap with the coordinate
    """
    imlist = [findimages(xc,yc,dire,**kwargs) for dire in dirlist if dire is not []]
    if flatten:
        return reduce(lambda a,b: a+b, imlist, [])
    else:
        return imlist

def get_cutouts(xcoord,ycoord,xwidth,ywidth, coordsys='galactic',
        ignore_imagetypes=('_area','power','angle','weight','residual','smooth','model','mask','noise','label'),
        flist='find', savedir=None, clobber=True):
    """
    Create cutouts from all possible images in the searched directories.
    """

    if flist == 'find':
        # Get list
        flist = find_all_images(xcoord, ycoord, dirs, coordsys=coordsys)
    elif type(flist) not in (tuple,list):
        raise TypeError("flist must be a list (or tuple) of filenames")

    # filter out image types I don't want
    for ign in ignore_imagetypes:
        flist = agpy.grep.grepv(ign, flist)

    # make all the cutouts
    cutouts = []
    for fn in flist:
        try:
            co = agpy.cutout.cutout(fn, xcoord, ycoord, xwidth, ywidth, units='wcs', coordsys=coordsys)
        except agpy.cutout.DimensionError:
            header = pyfits.getheader(fn)
            wcs = pywcs.WCS( agpy.cubes.flatten_header(header) )
            xc,yc,xw,yw = agpy.cubes.aper_world2pix((xcoord,ycoord,xwidth,ywidth),wcs,coordsys=coordsys,wunit='degree')
            co = agpy.cubes.subcube( pyfits.getdata(fn), xcoord, xwidth, ycoord, ywidth, header=header, return_HDU=True, widthunits='pixels' )
        except Exception as ex:
            print "%s failed: " % fn, ex
            flist.remove(fn)
        cutouts.append(co)
        if savedir is not None:
            co.writeto("%s/%s" % (savedir, os.path.split(fn)[1].replace(".fits","_cutout.fits")),clobber=clobber)

    return cutouts

def find_directories(rootdir):
    """
    Find directories containing FITS files
    """

    paths = []
    for root,dirs,files in os.walk(str(rootdir)):
        if ".fits" in (str(os.path.splitext(f)[1]) for f in files): 
            paths.append(root)

    return paths

default_rootdirs = ['/Volumes/disk2/data/','/Volumes/WD_2/HiGal/','/Volumes/disk4/']

dirs = reduce(lambda a,b: a+b, [find_directories(dr) for dr in default_rootdirs], []) 
