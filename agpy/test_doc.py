from pylab import *
import pylab
import agpy.test_doc

print "beta.__module__:",beta.__module__

for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'

