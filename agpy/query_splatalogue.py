import atpy
import urllib2
import tempfile

"""
Module to search Splatalogue.net via splat, modeled loosely on
ftp://ftp.cv.nrao.edu/NRAO-staff/bkent/slap/idl/
"""

length_dict = {'meters':1.0,'m':1.0,
        'centimeters':1e-2,'cm':1e-2,
        'millimeters':1e-3,'mm':1e-3,
        'nanometers':1e-9,'nm':1e-9,
        'micrometers':1e-6,'micron':1e-6,'microns':1e-6,'um':1e-6,
        'kilometers':1e3,'km':1e3,
        'angstroms':1e-10,'A':1e-10,
        }

def query_splatalogue(minwav,maxwav,waveunits='m',root_url='http://find.nrao.edu/splata-slap/slap'):
    """
    Acquire an atpy table of a splatalogue searched based on wavelength.

    Future work will allow queries based on other parameters.
    """

    if waveunits in length_dict:
        minwav = minwav * length_dict[waveunits]
        maxwav = maxwav * length_dict[waveunits]

    query_url = "%s?REQUEST=queryData&WAVELENGTH=%f/%f" % (root_url,minwav,maxwav)
    #query_url = urllib2.Request(url=root_url,
    #        data="queryData",{"WAVELENGTH":"%f/%f" % (minwav,maxwav)})

    U = urllib2.urlopen(query_url)

    # Normally it would be possible to do this:
    # t = atpy.Table(U,type='vo')
    # instead we have to write to a file and flush it. 

    R = U.read()
    U.close()
    tf = tempfile.NamedTemporaryFile()
    #for line in R:
    #    print >>tf,line.strip()
    print >>tf,R
    print tf.name
    tf.file.flush()
    t = atpy.Table(tf.name,type='vo')

    return t

