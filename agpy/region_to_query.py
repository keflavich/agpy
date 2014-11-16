#!/usr/bin/env python
"""
SIMBAD query a region
"""
import os
import sys
import ds9
import pyregion
import webbrowser
import numpy as np


xpa = sys.argv[1]
dd = ds9.ds9(xpa)

if len(sys.argv) > 2:
    regionid = int(sys.argv[2])
else:
    # use the most recent region
    regionid = -1

rstringlist = dd.get('regions -system wcs -prop select 1')
regions = pyregion.parse(rstringlist)
if len(regions) == 0:
    sys.exit("No regions found")


reg = regions[regionid]

coord_fmt = 'Gal' if reg.coord_format == 'galactic' else 'fk5'
coord = reg.coord_list[:2]

url = "http://simbad.u-strasbg.fr/simbad/sim-coo?CooFrame={frame}&Coord={coord}".format(frame=coord_fmt,
                                                                                        coord='{0} {1}'.format(*coord))

if reg.name == 'circle':
    radius = reg.coord_list[2]
    url += "&Radius={radius}&Radius.unit=deg".format(radius=radius)
else:
    radius = 2/60.

print url

webbrowser.open(url)

#dd.set('catalog simbad')
#dd.set('catalog coordinate {coord0} {coord1} {fmt} size {radius} {radius}'.format(coord0 = coord[0],
#                                                                                         coord1 = coord[1],
#                                                                                         fmt = coord_fmt,
#                                                                                         radius = radius))

from astroquery.simbad import Simbad
from astropy import units as u
from astropy import coordinates
from astropy.io import votable
import tempfile

C = coordinates.SkyCoord(coord[0]*u.deg, coord[1]*u.deg, frame=reg.coord_format)

tbl = Simbad.query_region(coordinates=C, radius=radius*u.deg)
tbl.pprint()

tbl['RA'] = [x.replace(" ",":") for x in tbl['RA']]
tbl['DEC'] = [x.replace(" ",":") for x in tbl['DEC']]
col = tbl['MAIN_ID']
tbl.remove_column('MAIN_ID')
tbl.add_column(col, index=2)
#tbl['MAIN_ID'] = [x.replace(" ","_").replace("[","").replace("]","") for x in tbl['MAIN_ID']]
bad = np.array([(x['RA'].count(":") < 2) or (x['DEC'].count(":") < 2) for x in tbl])
tbl = tbl[~bad]
tbl = tbl.filled(0)

for row_id,row in enumerate(tbl):
    for key,val in zip(row.dtype.names, row):
        if val == '':
            print row_id,key
            tbl[row_id][key] = 0

temp = tempfile.NamedTemporaryFile(suffix='.tsv')
tbl.write(temp.name, format='ascii.csv', delimiter='\t')
dd.set('catalog import tsv {0}'.format(temp.name))
dd.set('catalog sky fk5')
dd.set('catalog psky fk5')
#except:
#    vot = votable.table.from_table(tbl)
#    temp = tempfile.NamedTemporaryFile(suffix='.xml')
#    vot.to_xml(temp.name)
#    dd.set('catalog load {0}'.format(temp.name))
#
