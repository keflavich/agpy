"""
Make cutouts with all available data of a given position
"""
import glob
import pyfits
import pywcs
import coords
import os
import agpy.cubes

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
    wcs = pywcs.WCS(header)

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

def get_cutouts(xcoord,ycoord,coordsys='galactic'):

    pass

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
