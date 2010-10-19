from agpy import cubes
import numpy as np

zz,yy,xx = np.indices([50,50,50])
rr = np.sqrt(25-((zz-25)**2+(yy-25)**2+(xx-25)**2))

smcube = cubes.smooth_cube(rr,ignore_nan=True)

print "If you've gotten here, parallel smoothing works, but there's a lot of testing left to be done"
