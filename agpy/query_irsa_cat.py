import urllib2
import time
import os
import sys
import atpy

def query_irsa_cat( targetname_OR_coords, catalog_name='fp_psc', radius=60,
        radunits='arcsec', outfile=None, debug=False, overwrite=False,
        outfmt='ipac'):
    """
    #;+
    #; NAME: 
    #;    QUERY_IRSA_CAT
    #;
    #; PURPOSE: 
    #;    Query any catalog in the NASA/IPAC Infrared Science Archive (IRSA)
    #;    database by position or resolvable target name.
    #; 
    #; EXPLANATION:
    #;    Uses the IDL SOCKET command to provide a query of any catalog 
    #;    in the IRSA (http://irsa.ipac.caltech.edu/) database over the Web and
    #;    return results in an IDL structure.  If outfile is set, it saves
    #;    the query as an IPAC table.  This can be slow for large query
    #;    results, so only write a file if needed.    
    #;     
    #; CALLING SEQUENCE: 
    #;    info = query_irsa_cat(targetname_or_coords, [catalog=catalog, radius=radius,
    #;                           radunits=radunits, outfile=outfile])
    #;
    #; INPUTS: 
    #;
    #;    TARGETNAME_OR_COORDS - Either a string giving a resolvable target
    #;           name (with J2000 coordinates determined by NED or SIMBAD), or a 
    #;           2-element numeric vector giving the J2000 right ascension 
    #;           and declination, both in degrees.
    #;
    #; OPTIONAL INPUT:
    #;
    #;    CATALOG - string giving the identifier string of the IRSA catalog to be
    #;           searched.  The complete list of catalogs and identifier strings is available in
    #;           XML format at:
    #;             http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mode=xml
    #;           or as an IPAC Table (ascii) at:
    #;             http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-scan?mode=ascii
    #;
    #;           In the table, the identifier string is in the "catname" column.
    #;
    #;           If this keyword is omitted, the program will query the 2MASS point
    #;           source catalog
    #;
    #;           Examples of IRSA catalogs include:  
    #;              'fp_psc' - 2MASS Point Source Catalog
    #;              'iraspsc' - IRAS Point Source Catalog v2.1 (PSC)
    #;              'irasfsc' - IRAS Faint Source Catalog v2.0
    #;              'cosmos_ib_phot' - COSMOS Intermediate and Broad Band Photometry Catalog 2008
    #;              'akari_irc' - Akari/IRC Point Source Catalogue
    #;
    #;    RADIUS - scalar input of the radius of the search.  By default it has a value of 60. IRSA
    #;           catalogs have maximum allowable search radii.  These are listed on the corresponding
    #;           web interface page for the catalog search, or in the nph-scan return table in the
    #;           "coneradius" column.
    #;    
    #;    RADUNITS - string giving the units of the radius.  By default it is 'arcsec'.
    #;    
    #;    OUTFILE - if present, the search results are written to a file with this name.
    #;
    #;    DEBUG - /DEBUG provides some additional output.
    #;
    #; OUTPUTS: 
    #;    info - Anonymous IDL structure containing information on the catalog.  The structure
    #;           tag names are taken from the catalog column names.  If no objects were found in 
    #;           the catalog, the structure values will be empty or zero.  If any input parameters
    #;           (e.g. catalog name) is invalid, the structure will have no
    #;           content fields other than info.CREATED.
    #;
    #;           If the query fails or is invalid, the function returns a value of -1.  
    #;
    #; EXAMPLES: 
    #;           (1) Plot a histogram of the J magnitudes of all 2MASS
    #;               point sources within 10 arcminutes of the center of the
    #;               globular cluster M13.  Save the IPAC table. 
    #;
    #;             IDL> info = query_irsa_cat('m13',radius=10,radunits='arcmin',outfile='save.tbl')
    #;             IDL> help,/struct,info
    #;             IDL> plothist,info.j_m,xran=[10,20]
    #;
    #;           (2) Find the position of the faintest IRAS 60 micron
    #;               source within 1 degree of central position of the
    #;               COSMOS survey (10h00m28.6s +02d12m21.0s in J2000)
    #;  
    #;             IDL> info = query_irsa_cat([150.11917,2.205833], catalog='irasfsc', radius=1, radunits='deg')
    #;             IDL> help,/struct,info
    #;             IDL> idx = where(info.fnu_60 eq min(info.fnu_60))
    #;             IDL> print, (info.ra)[idx], (info.dec)[idx]
    #;
    #; PROCEDURES USED:
    #;    READ_IPAC_VAR  available from IRSA
    #;    IS_IT_NUMBER   available from IRSA
    #;    WEBGET()       from IDLastro.  Must be 2009 update or later.
    #;
    #; NOTES:
    #;    The program writes an output IPAC table file only if the
    #;    OUTFILE keyword is set.
    #;
    #; MODIFICATION HISTORY:
    #;    Harry Teplitz   September 2010 (based closely on queryvizier.pro)
    #;    Removed requirement of writing/reading IPAC table file - TYB, May 2011
    """

    if outfile is not None:
        writefile=outfile
        check = os.path.exists(writefile)
        if check and not overwrite:
            raise IOError("ERROR: File %s exists and overwrite=False." % (writefile))

    ##;;;;;;;;;;;;;;;;;;;  CONSTRUCT THE PARTS OF THE QUERY STRING

    root = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query'

    #;;;; CATALOG STRING

    catstr='&catalog='+catalog_name

    #;;;; OBJECT STRING

    target = targetname_OR_coords

    if len(target) == 2:
        ra = float(target[0])
        dec = float(target[1])
        #;  IF dec GT 0 THEN sign='+' ELSE sign=''
        objstr = '&objstr=%f+%f' % (ra,dec)
    else:
        object = target.replace("+","%2B")
        object = object.strip().replace(' ','+')
        objstr = '&objstr='+object

        if len(objstr) <= 8:
            raise ValueError("ERROR: Empty string not allowed for target.")

    #;;;;  SEARCH SHAPE AND SIZE

    spatial_str='Cone' 
    spatial_param_name=['radius','radunits']
    spatial_param_value_str = [str(radius), radunits]

    nspat = len(spatial_param_name)

    spatstr = '&spatial='+spatial_str
    spatparstr = ''

    for ii in xrange(nspat):
        spatparstr += '&'+spatial_param_name[ii]+'='+spatial_param_value_str[ii]
            
    # specify output table type; default is IPAC
    outfmt_dict = {'ipac':1,'vo':3}
    if outfmt not in outfmt_dict:
        raise ValueError("Output format %s not recognized." % outfmt)

    out_fmt = '?outfmt=%i' % outfmt_dict[outfmt]

    #;;;; combine into query string

    url_q = root+out_fmt+objstr+spatstr+spatparstr+catstr
    if debug:
        print url_q

    # just kept for debug purposes
    # url_q_old = 'http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?outfmt=3&objstr=10.68469-41.26904&spatial=Cone&radius=100&radunits=arcsec&catalog=fp_psc'

    if debug:
        t0 = time.time()

    # get the URL 
    request = urllib2.Request(url_q)
    ret = urllib2.urlopen(request)
    if ret.msg != 'OK':
        print "Message: "+ret.msg

    if debug:
        print time.time()-t0

    # use ATPY to parse the table
    table = atpy.Table(ret, type=outfmt)

    #;;;;;  If requested, write the output to the outputfile
    if outfile is not None:
        table.writeto(writefile,type=outfmt)

    return table
