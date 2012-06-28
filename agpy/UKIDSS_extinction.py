from astrodata import ukidss,magpis
import matplotlib
import matplotlib._cntr as _cntr
import pylab
import numpy as np
try:
    import astropy.wcs as pywcs
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
    import pywcs

def get_data(glon,glat,radius=20,save=True,overwrite=False,get_images=True):
    """
    Download UKIDSS data at specified GLON/GLAT with specified size

    Parameters
    ----------
    glon : float
    glat : float
        Galactic latitude and longitude at the center
    radius : float
        Size of cutout (symmetric) in arcminutes
    save : bool
        Save FITS and catalog files?
    get_images : bool
        Download the images in addition to the catalog?
    """

    R = ukidss.UKIDSSQuery()
    if get_images:
        images = R.get_image_gal(glon,glat,size=radius,save=save,overwrite=overwrite)
    cat = R.get_catalog_gal(glon,glat,radius=radius,save=save,overwrite=overwrite)[0][1]
    cleancat = ukidss.ukidss.clean_catalog(cat)

    return cleancat


def make_densitymap(cat, pixsize=7.2, save_prefix="densmap_", overwrite=False):
    """
    Create point source density maps in glon/glat

    Parameters
    ----------
    cat : astropy.io.fits.hdu.table.BinTableHDU
        The catalog to densimap
    pixsize : float 
        Pixel size in which to bin data (arcseconds)
    save_prefix : string
        Prepends to the .fits filename that will be saved; band name will be
        appended
    """

    jhk = []
    for band in ('J','H','K_1'):
        mask = ((cat[band+'CLASS']!=-9999) * (cat[band+'ERRBITS'] <
            42) * (cat[band+'ERRBITS'] > -1) * ((cat['PRIORSEC'] ==
                cat['FRAMESETID']) + (cat['PRIORSEC']==0))
            * (cat[band+'PPERRBITS']!=64) * (cat[band+'PPERRBITS'] < 60) #70000)
            )
        lonK=cat['L'][mask]
        latK=cat['B'][mask]

        binsize = pixsize
        mapsize = np.array( [(lonK.max()-lonK.min()),(latK.max()-latK.min())] )
        mapsize_pix = (mapsize)/(binsize/3600.)
        H,histlon,histlat = np.histogram2d(lonK,latK,bins=mapsize_pix)

        F = pyfits.PrimaryHDU()
        F.header.update('CRPIX1',mapsize_pix[1]/2.+1)
        F.header.update('CRPIX2',mapsize_pix[0]/2.+1)
        F.header.update('CRVAL1',np.median(histlon) )
        F.header.update('CRVAL2',np.median(histlat) )
        F.header.update('CDELT1', binsize/3600.     )
        F.header.update('CDELT2', binsize/3600.     )
        F.header.update('CTYPE1', 'GLON-CAR'     )
        F.header.update('CTYPE2', 'GLAT-CAR'     )
        F.data = H.T
        F.writeto('%s%s.fits' % (save_prefix,band),clobber=overwrite)
        jhk.append(F)

    return jhk

def get_image(glon,glat,radius=20,save=True,overwrite=False):
    """
    Download Bolocam image data at specified GLON/GLAT with specified size

    Parameters
    ----------
    glon : float
    glat : float
        Galactic latitude and longitude at the center
    radius : float
        Size of cutout (symmetric) in arcminutes
    save : bool
        Save FITS?
    """

    return magpis.get_magpis_image_gal(glon, glat, size=radius, save=save, overwrite=overwrite)

def get_contours(fitsfile, av=10.):
    """
    Given a Bolocam FITS file, return the contours at a given flux level
    """

    header = fitsfile[0].header
    img = fitsfile[0].data

    # from Foster 2012
    av_to_jy = 6.77e22/9.4e20 # cm^-2 / Jy / (cm^-2 / AV) = AV/Jy
    #if header.get('BGPSVERS').strip()=='1.0':
    av_to_jy /= 1.5

    contour_level = av / av_to_jy

    wcs = pywcs.WCS(header)
    #wcsgrid = wcs.wcs_pix2world( np.array(zip(np.arange(wcs.naxis1),np.arange(wcs.naxis2))), 0 ).T
    yy,xx = np.indices(img.shape)

    img[img!=img] = 0
    C = _cntr.Cntr(yy,xx,img)
    paths = [p for p in C.trace(contour_level) if p.ndim==2]

    wcs_paths = [wcs.wcs_pix2world(p,0) for p in paths]

    return wcs_paths

def histeq(im,nbr_bins=256):

   #get image histogram
   imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = 255 * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   im2 = np.interp(im.flatten(),bins[:-1],cdf)

   return im2.reshape(im.shape)#, cdf


def show_contours_on_extinction(contours, jhk, color='k'):
    """
    Given contours from get_contours and a list of JHK images from make_densitymap, plot things
    """

    header = jhk[0].header
    rgb = ([im.data for im in jhk])
    alpha = np.array(rgb).sum(axis=0)
    alpha /= alpha.max()
    alpha *= 0.5
    alpha += 0.5
    alpha[:]=1
    #alpha = histeq(alpha)
    rgb.append(alpha)
    rgb = np.array(rgb).T
    rgb[:,:,:3] /= 10.

    wcs = pywcs.WCS(header)
    xglon,yglat = wcs.wcs_pix2world( np.array(zip(np.arange(wcs.naxis1),np.arange(wcs.naxis2))), 0 ).T

    pylab.imshow(rgb,extent=[xglon.min(),xglon.max(),yglat.min(),yglat.max()])
    for C in contours:
        pylab.plot(*C.T.tolist(),color=color)

def contour_segments(p):
    return zip(p, p[1:] + [p[0]])

def contour_area(p):
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(p)))


if __name__ == "__main__":

    if not 'glon' in locals():
        glon = 10.62
        glat = -0.38
        radius = 20
    cat = get_data(glon,glat,radius=radius,overwrite=True,get_images=False)
    jhk = make_densitymap(cat,overwrite=True,save_prefix="G%07.3f%+08.3f_densmap_" % (glon,glat))
    bgps = get_image(glon,glat,radius=radius,overwrite=True)
    contours = get_contours(bgps)
    show_contours_on_extinction(contours, jhk)
