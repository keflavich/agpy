#!python

import math
import warnings
warnings.filterwarnings('ignore','masked')
import pylab
from pylab import *
for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'
import matplotlib
import pyfits
import numpy 
from mad import MAD,nanmedian
#from matplotlib import patches
from matplotlib.patches import Rectangle,FancyArrow,Circle,Ellipse
from matplotlib.lines import Line2D
from matplotlib.widgets import Cursor, MultiCursor
import matplotlib.cm as cm
#from Scientific.IO import NetCDF
#from scipy.io import netcdf
import time
import re
import os
import subprocess
import copy
import idlsave
import gaussfitter
import mpfit
from PCA_tools import *
from AG_image_tools.drizzle import drizzle
from agpy import smooth
from guppy import hpy
heapy = hpy()

matplotlib.rcParams['image.origin']='lower'
matplotlib.rcParams['image.interpolation']='nearest'
matplotlib.rcParams['image.aspect']=1
matplotlib.rcParams['axes.color_cycle'] = [list(clr) for clr in matplotlib.cm.brg(linspace(0,1,144))]
matplotlib.defaultParams['image.origin']='lower'
matplotlib.defaultParams['image.interpolation']='nearest'
matplotlib.defaultParams['image.aspect']=1
if matplotlib.rcParams['text.usetex']: texOn = True
else: texOn = False
# matplotlib.rcParams['text.usetex']=False
# matplotlib.defaultParams['text.usetex']=False

class lazydata(object):
    def __init__(self, varname, structname='bgps', reshape=None, flag=True):
        self.varname = varname
        self.structname = structname
        self.reshape = reshape
        self.flag = flag

    def __getattr__(self,attribute):
        return getattr(self,attribute)

    def __get__(self, obj, type=None):
        t0 = time.time()
        if obj.__dict__.has_key(self.varname):
            # technically I don't think this should ever be called
            # but it is. ARGH.
            #print "Getting (instead of computing) %s" % self.varname
            return obj.__dict__[self.varname]
        else:
            print "Computing %s " % self.varname
            if self.flag:
                obj.__dict__[self.varname] = obj.__dict__[self.structname][self.varname][0][obj.whscan,:].astype('float')
                obj.__dict__[self.varname][obj.whempty,:] = NaN
                obj.__dict__[self.varname].shape = obj.datashape
                obj.__dict__[self.varname] = nantomask(obj.__dict__[self.varname])
                try:
                    obj.__dict__[self.varname].mask[obj.flags > 0] = True
                except TypeError:
                    obj.__dict__[self.varname].mask = (obj.flags > 0)
            else:
                # flagging is extra work, skip it but still do the reshaping aspects
                obj.__dict__[self.varname] = obj.__dict__[self.structname][self.varname][0][obj.whscan,:].astype('float')
                obj.__dict__[self.varname].shape = obj.datashape

            # not used if self.reshape is not None:
            # not used     obj.__dict__[self.varname] = reshape( obj.mapstr[self.varname][0][obj.whscan,:] , obj.datashape )
            print "Finished computing %s in %0.3g seconds" % (self.varname,time.time()-t0)
            return obj.__dict__[self.varname]

class Flagger:
    """
    Write out a file with appropriate flagging commands for use in IDL / later editing
    Example:

        import pyflagger
        f = pyflagger.Flagger('050906_o11_raw_ds5.nc_indiv13pca_timestream00.fits','050906_o11_raw_ds5.nc')
        f.plotscan(0)
        f.close()

    Key commands:
      left click - flag
      right click - unflag
      n - next scan
      p - previous scan
      q - save and quit
      Q - quit (no save)
      . - point to this point in the map
      f - plot footprint of array at this time point
      R - reverse order of flag boxes (to delete things hiding on the bottom)
      r - redraw
      d - delete flag box
      t - flag timepoint
      s - flag scan
      w - flag Whole scan (this is the same as s, except some python backends catch / steal 's')
      S - unflag scan
      b - flag bolometer
      T - unflag timepoint
      B - unflag bolometer
      c - toggle current scan
      v - display data value
      P - display the PCA decomposition of the displayed timestream
      o - make a map of the array at the sampled time
      z - display the power spectra of the displayed timestream (use 'C' to plot one)
      Z - display the power spectra of the displayed timestream over all time
      C,L - plot Column/Line
      j - plot whole timestream for selected bolo
      a - create a footprint movie between two selected points
      M,m - flag highest, lowest point in map
 
    Map Key Commands:
      c - toggle current scan
      . - show point in timestream
      click - show point in timestream
      middle click - list all points that contribute to that pixel
      r - redraw

    """



    def __init__(self, filename, debug=False, npca=13, **kwargs):
        # Initialize plots first
        pylab.figure(0)
        pylab.figure(1,figsize=[16,12])
        pylab.figure(2,figsize=[16,12])
        self.filename = filename
        self.debug = debug
        self.npca = npca
        if filename[-4:] == 'fits':
            self._loadfits(filename,**kwargs)
        elif filename[-3:] == 'sav':
            self._loadsav(filename,**kwargs)

        self.help = """

    Key commands:
      left click - flag
      right click - unflag
      n - next scan
      p,N - previous scan
      q - save and quit
      Q - quit (no save)
      . - point to this point in the map
      f - plot footprint of array at this time point
      R - reverse order of flag boxes (to delete things hiding on the bottom)
      r - redraw
      d - delete flag box
      t - flag timepoint
      s - flag scan
      w - flag Whole scan (this is the same as s, except some python backends catch / steal 's')
      S,W - unflag scan
      b - flag bolometer
      T - unflag timepoint
      B - unflag bolometer
      c - toggle current scan
      v - display data value
      P - display the PCA decomposition of the displayed timestream
      o - make a map of the array at the sampled time
      z - display the power spectra of the displayed timestream (use 'C' to plot one)
      Z - display the power spectra of the displayed timestream over all time
      C,L - plot Column/Line
      j - plot whole timestream for selected bolo
      a - create a footprint movie between two selected points
      M,m - flag highest, lowest point in map
      e - expsub current plane
 
    Map Key Commands:
      c - toggle current scan
      . - show point in timestream
      click - show point in timestream
      middle click - list all points that contribute to that pixel
      r - redraw

        """

    def _loadfits(self, filename, ncfilename='', flagfile='', mapnum='', axis=None, **kwargs):
        fnsearch = re.compile(
                '([0-9]{6}_o[0-9b][0-9]_raw_ds5.nc)(_indiv[0-9]{1,2}pca)').search(filename)
        ncsearch = re.compile(
                '[0-9]{6}_o[0-9b][0-9]_raw_ds5.nc').search(ncfilename)
        if fnsearch is None:
            print "Couldn't find the correct prefix in the filename" \
                    +" - expected form like 050906_o11_raw_ds5.nc_indiv13pca_timestream00.fits"
            return
        mapnumsearch = re.compile('([0-9]{2})(\.fits)').search(filename)
        if mapnumsearch is not None and mapnum=='':
            mapnum = mapnumsearch.groups()[0]
        else:
            mapnum = '01'

        if fnsearch.groups()[0] == ncsearch.group():
            self.ncfilename = ncfilename # self.pathprefix+fnsearch.groups()[0]
            self.readncfile()
        else:
            print "Warning: the NCDF filename doesn't match the input fits file name."\
                    + "You'll probably get errors and your work won't be saved."
            self.ncfilename = self.pathprefix+fnsearch.groups()[0]

        self.fileprefix = fnsearch.group()
        self.pathprefix = filename[:fnsearch.start()]
        self.tsfn = self.pathprefix+self.fileprefix+"_timestream00.fits"
        self.tsfile = pyfits.open(self.tsfn)
        self.mapfn = self.pathprefix+self.fileprefix+"_map"+mapnum+".fits"
        self.mapfile = pyfits.open(self.mapfn)
        self.map = self.mapfile[0].data
        self.map[numpy.isnan(self.map)] = 0
        self.tstomapfn = self.pathprefix+self.fileprefix+"_tstomap.fits"
        self.tstomapfile = pyfits.open(self.tstomapfn)
        self.tstomap = self.tstomapfile[0].data

#      self.outfile = open(self.pathprefix+"log_"+self.fileprefix+"_flags.log",'a')
        self.data = self.tsfile[0].data
        self.flagfn = self.pathprefix+self.fileprefix+"_flags.fits"
#      if os.path.exists(self.flagfn):
#          self.flagfile = pyfits.open(self.flagfn)
#          self.flags = self.flagfile[0].data
#          if self.flags.shape != self.data.shape:
#              print "Flags / data shape are different.",self.flags.shape,self.data.shape
#      else:
#          self.flagfile = copy.copy(self.tsfile)
#          self.flags = zeros(self.data.shape,dtype='int')
        self._initialize_vars(**kwargs)
    
    def _initialize_vars(self,vmax=None):
        #print >>self.outfile,"Started a new session at "+time.asctime()
        self.reset()
        self.counter = 0
        self.mouse_up = False
        self.connected = 0
        #self.renderer = matplotlib.backends.backend_agg.RendererAgg
        print "There are %i scans" % (self.data.shape[0])
        self.maxscan = self.data.shape[0]
        self.rectangles=[[] for i in xrange(self.maxscan)]
        self.lines=[[] for i in xrange(self.maxscan)]
        self.arrows=[]
        self.maparrows=[]
        self.connections=[]
        self.mapconnections=[]
        self.md = 0
        self.mu = 0
        self.key = 0
        self._lastkey = None
        self.scannum = 0
        self.fignum = 1
        self.open = 1
        self.currentscan = 0
        self.aspect = float(self.data.shape[2])/float(self.data.shape[1])
        self.plotfig = None
        self.bolofig = None
        self.mapfig = None
        self.flagfig = None
        self.datafig = None
        self.scanim = None
        self.PCAflag = False
        self.powerspec_plotted = False
        self.powerspectra_whole = None
        self.gaussfit=None

        self.showmap(vmax=vmax)
        self.dcon()
    
    def _loadsav(self, savfile, flag=True, **kwargs):
        memtot = heapy.heap().size / 1024.0**3
        print "Beginning IDLsave file read. %0.3g GB used" % memtot
        t0 = time.time()
        sav = idlsave.read(savfile)
        memtot = heapy.heap().size / 1024.0**3
        t1 = time.time()
        print "Finished reading IDLsave file in %i seconds using %0.3g GB" % (t1 - t0,memtot)
        self.bgps = sav.get('bgps')
        memtot = heapy.heap().size / 1024.0**3
        t2 = time.time()
        print "Set bgps variable in %i seconds using %0.3g GB" % (t2 - t1,memtot)
        self.mapstr = sav.get('mapstr')
        self.needed_once_struct = sav.get('needed_once_struct')
        if self.needed_once_struct is None:
            neededoncefile = savfile.replace('preiter','neededonce').replace('postiter','neededonce')
            if os.path.exists(neededoncefile):
                sav_once = idlsave.read(neededoncefile)
                self.needed_once_struct = sav_once.get('needed_once_struct')
        t3 = time.time()
        print "Completed IDLsave file read in %f seconds." % (t3 - t0)

        self.ncfilename = savfile
        self.tsfile = None

        fnsearch = re.compile(
                '([0-9]{6}_o[0-9b][0-9]_raw_ds[125].nc)(_indiv[0-9]{1,2}pca)').search(savfile)
        self.fileprefix = fnsearch.group()
        self.pathprefix = savfile[:fnsearch.start()]

        self.ncscans = self.bgps['scans_info'][0]
        self.sample_interval = self.bgps['sample_interval'][0]
        if len(self.ncscans.shape) == 1: self.ncscans.shape = [1,2]
        self.scanlengths = self.ncscans[:,1]+1-self.ncscans[:,0]
        self.scanlen = numpy.max(self.scanlengths)
        self.ncflags = self.bgps['flags'][0] 
        self.timelen = self.ncflags.shape[0]
        self.nbolos = self.ncflags.shape[1]
        self.nscans = self.ncscans.shape[0]
        self.ncbolo_params = self.bgps['bolo_params'][0]
        self.ncbolo_indices = self.bgps['bolo_indices'][0]
        #self.bolo_indices = asarray(nonzero(self.ncbolo_params[:,0].ravel())).ravel()
        self.bolo_indices = self.ncbolo_indices
        self.ngoodbolos = self.bolo_indices.shape[0]
        self.whscan = asarray([arange(self.scanlen)+i for i,j in self.ncscans[:,:2]]).ravel()
        self.scanstarts = arange(self.nscans)*self.scanlen
        self.whempty = concatenate([arange(i+j,i+self.scanlen) for i,j in zip(self.scanstarts,self.scanlengths) ]).ravel()
        self.whscan[self.whempty] = 0

        self.tsshape = [self.nscans*self.scanlen,self.ngoodbolos]
        self.datashape = [self.nscans,self.scanlen,self.ngoodbolos]

        t4 = time.time()
        memtot = heapy.heap().size / 1024.0**3
        print "Beginning array reshaping with %f seconds elapsed, %0.3g GB used." % (t4 - t0, memtot)

        #class lazy_whscan(object):
        #    class data(object):
        #        def __get__(self, obj, varname, type=None):
        #            print "Computing ... "
        #            obj.__dict__[varname] = self.bgps[varname][0][self.whscan,:].astype('float')
        #    data = data()

        #self.flags          = lazydata('flags',flag=False) # self.bgps['flags'][0][self.whscan,:]
        setattr(self.__class__, 'flags', lazydata('flags',flag=False))
        self.flags.shape    = self.datashape

        if self.needed_once_struct is not None:
            print "Loading 'raw' and 'dc_bolos' from needed_once_struct"
            #self.raw         = lazydata('raw', 'needed_once_struct') #self.needed_once_struct['raw'][0][self.whscan,:].astype('float')
            #self.dc_bolos    = lazydata('dc_bolos', 'needed_once_struct') #self.needed_once_struct['dc_bolos'][0][self.whscan,:].astype('float')
            setattr(self.__class__, 'raw', lazydata('raw', 'needed_once_struct',flag=flag)) #self.needed_once_struct['raw'][0][self.whscan,:].astype('float')
            setattr(self.__class__, 'dc_bolos', lazydata('dc_bolos', 'needed_once_struct',flag=flag)) #self.needed_once_struct['dc_bolos'][0][self.whscan,:].astype('float')
        elif self.bgps.dtype.fields.has_key('raw'):
            print "Loading 'raw' and 'dc_bolos' from bgps"
            #self.raw         = lazydata('raw') # self.bgps['raw'][0][self.whscan,:].astype('float')
            #self.dc_bolos    = lazydata('dcbolos') # self.bgps['dc_bolos'][0][self.whscan,:].astype('float')
            setattr(self.__class__, 'raw', lazydata('raw',flag=flag)) #self.needed_once_struct['raw'][0][self.whscan,:].astype('float')
            setattr(self.__class__, 'dc_bolos', lazydata('dc_bolos',flag=flag)) #self.needed_once_struct['dc_bolos'][0][self.whscan,:].astype('float')
        #self.astrosignal    = self.bgps['astrosignal'][0][self.whscan,:].astype('float')
        #self.atmosphere     = self.bgps['atmosphere'][0][self.whscan,:].astype('float')
        #self.ac_bolos       = self.bgps['ac_bolos'][0][self.whscan,:].astype('float')
        #self.atmo_one       = self.bgps['atmo_one'][0][self.whscan,:].astype('float')
        #self.noise          = self.bgps['noise'][0][self.whscan,:].astype('float')
        #self.scalearr       = self.bgps['scalearr'][0][self.whscan,:].astype('float')
        self.scale_coeffs   = self.bgps['scale_coeffs'][0].astype('float')
        #self.weight         = self.bgps['weight'][0][self.whscan,:].astype('float')
        #self.zeromedian     = self.astrosignal * 0
        t5 = time.time()
        memtot = heapy.heap().size / 1024.0**3
        print "Finished array reshaping in %f seconds, %0.3g GB used." % (t5 - t4, memtot)
        print "Beginning array flagging."

        #try:
        #    self.mapped_astrosignal = self.bgps['mapped_astrosignal'][0][self.whscan,:].astype('float')
        #except ValueError:
        #    self.mapped_astrosignal = copy.copy(self.astrosignal)

        datums=['astrosignal','atmosphere','ac_bolos','atmo_one','noise','scalearr','weight','mapped_astrosignal']
        for d in datums:
            if hasattr(self.bgps,d): # for version one, may not have some...
                setattr(self.__class__, d, lazydata(d,flag=flag))
            #self.__dict__[d] = lazydata(d)
            #if self.__dict__.has_key(d):
            #    self.__dict__[d][self.whempty,:] = NaN
            #    self.__dict__[d].shape = self.datashape
            #    self.__dict__[d] = nantomask(self.__dict__[d])
            #    try:
            #        self.__dict__[d].mask[self.flags > 0] = True
            #    except TypeError:
            #        self.__dict__[d].mask = (self.flags > 0)
        #import pdb; pdb.set_trace()
        self.weight_by_bolo = self.weight.mean(axis=0).mean(axis=0)
        if hasattr(self.bgps,'mapped_astrosignal'):
            setattr(self.__class__, 'mapped_astrosignal', lazydata('mapped_astrosignal',flag=flag))
        else:
            setattr(self.__class__, 'mapped_astrosignal', lazydata('astrosignal',flag=flag))

        if list(self.ac_bolos.shape) != self.datashape:
            import pdb; pdb.set_trace()
        self.data = self.ac_bolos

        self.ncfile = None
        self.flagfn = savfile.replace("sav","_flags.fits")

        self.map      = nantomask( self.mapstr['astromap'][0] )
        self.default_map = nantomask( self.mapstr['astromap'][0] )
        self.model    = nantomask( self.mapstr['model'][0] )
        self.noisemap = nantomask( self.mapstr['noisemap'][0] )
        #setattr(self.__class__, 'tstomap', lazydata('ts',reshape=True, structname='mapstr'))
        #self.tstomap  = lazydata('ts',reshape=True) # 

        if not hasattr(self.bgps,'atmo_one'):
            print "Reading file as a v1.0.2 sav file"
            self.atmo_one = self.ac_bolos - self.astrosignal
            self.mapped_timestream = self.ac_bolos - self.atmosphere # apparently?
            self.scalearr = numpy.ones(self.datashape[1])[newaxis,:,newaxis]*self.scale_coeffs.swapaxes(0,1)[:,newaxis,:]
            self.version = 'v1.0'
        else:
            self.mapped_timestream = self.atmo_one - self.atmosphere + self.astrosignal
            self.version = 'v2.0'

        if self.map.sum() == 0:
            self.map  = nantomask( self.mapstr['rawmap'][0] )

        self.header = pyfits.Header(_hdr_string_list_to_cardlist( self.mapstr['hdr'][0] ))

        t6 = time.time()
        memtot = heapy.heap().size / 1024.0**3
        print "Finished array flagging in %f seconds, %0.3g GB used." % (t6 - t5, memtot)

        # don't delay this
        self.tstomap = reshape( self.mapstr['ts'][0][self.whscan,:] , self.datashape )
        t7 = time.time()
        memtot = heapy.heap().size / 1024.0**3
        print "Computed tstomap in %f seconds, %0.3g GB used." % (t7 - t6, memtot)

        self._initialize_vars(**kwargs)

        self.tsplot_dict = {'astrosignal': lambda: self.astrosignal if self.astrosignal.sum() != 0 else 0,
        'dc_bolos': lambda: self.dc_bolos*self.scalearr,
        'dc_bolos_noscale': lambda: self.dc_bolos,
        'dcbolos': lambda: self.dc_bolos*self.scalearr,
        'dcbolos_noscale': lambda: self.dc_bolos,
        'acbolos_noscale': lambda: self.ac_bolos,
        'ac_bolos_noscale': lambda: self.ac_bolos,
        'atmo_one': lambda:self.atmo_one,
        'acbolos': lambda:self.ac_bolos*self.scalearr,
        'ac_bolos': lambda:self.ac_bolos*self.scalearr,
        'atmosphere': lambda:self.atmosphere,
        'skysub_noscale': lambda:self.atmo_one - self.atmosphere,
        'new_astro': lambda: self.atmo_one - self.atmosphere,
        'new_astro_v1': lambda: self.lookup('PCA_astro_v1'),
        'residual': lambda:self.atmo_one - self.atmosphere - self.noise,
        'skysub': lambda:self.atmo_one - self.atmosphere + self.astrosignal,
        'default': lambda:self.atmo_one - self.atmosphere + self.astrosignal,
        'last_astrosignal': lambda:self.atmo_one - self.atmosphere - self.noise + self.astrosignal,
        'acbMatmo': lambda: self.ac_bolos - self.atmo_one - self.atmosphere,
        'acbMatmosphere': lambda: self.ac_bolos - self.atmosphere,
        'acbMatmoone': lambda: self.ac_bolos - self.atmo_one,
        'scale': lambda: self.scalearr,
        'weight': lambda: self.weight,
        'raw': lambda: self.raw,
        'rawscaled': lambda: self.raw * self.scalearr,
        'noise': lambda: self.noise,
        'newnoise': lambda: self.PCA_astro + self.astrosignal - self.astrosignal_from_model,
        'mapped_astrosignal': lambda: self.mapped_astrosignal,
        'mapped_timestream': lambda: self.mapped_timestream,
        'astrosignal_from_map': lambda: self.default_map.flat[self.tstomap],
        'astrosignal_from_model': lambda: self.model.flat[self.tstomap],
        'itermedian': lambda: itermedian(self.ac_bolos * self.scalearr),
        'zeromedian': lambda: self.atmo_one,
        'atmo_one_itermedian': lambda: itermedian(self.atmo_one),
        #'expsub': lambda: exponent_sub(self.lookup('atmo_one_itermedian')),
        'atmos_remainder': lambda: self.lookup('atmo_one_itermedian'),#self.lookup('expsub'),
        'atmos_remainder_v1': lambda: itermedian(self.atmo_one,scale=1.0,niter=1),
        'expmodel': lambda: self.lookup('atmo_one_itermedian') - self.lookup('expsub'),
        'first_sky': lambda: self.atmo_one - self.lookup('atmos_remainder'),
        'first_sky_v1': lambda: self.atmo_one - self.lookup('atmos_remainder_v1'),
        'astrosignal_premap': lambda: self.lookup('PCA_astro')+self.astrosignal,
        'PCA_atmo_v1':     lambda: reshape(unpca_subtract(numpy.nan_to_num(reshape(self.lookup('atmos_remainder_v1'),self.tsshape)),self.npca),self.datashape),
        'PCA_astro_v1':     lambda: reshape(pca_subtract(numpy.nan_to_num(reshape(self.lookup('atmos_remainder_v1'),self.tsshape)),self.npca),self.datashape),
        'PCA_atmo':     lambda: reshape(unpca_subtract(numpy.nan_to_num(reshape(self.lookup('atmos_remainder'),self.tsshape)),self.npca),self.datashape),
        'PCA_astro':     lambda: reshape(pca_subtract(numpy.nan_to_num(reshape(self.lookup('atmos_remainder'),self.tsshape)),self.npca),self.datashape),
        'PCA_astrosignal':   lambda: reshape(efuncs(reshape(self.astrosignal,self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_acb':     lambda: reshape(efuncs(reshape(self.ac_bolos,self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_zeromedian': lambda: reshape(efuncs(reshape(self.atmo_one,self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_itermedian': lambda: reshape(efuncs(reshape(self.lookup('itermedian'),self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_noise':   lambda: reshape(efuncs(reshape(self.noise,self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_default': lambda: reshape(efuncs(reshape(self.atmo_one - self.atmosphere + self.astrosignal,self.tsshape)),self.datashape) / self.nbolos**0.5,
        'PCA_atmos_remainder': lambda: reshape(efuncs(reshape(numpy.nan_to_num(self.lookup('atmos_remainder')),self.tsshape)),self.datashape) / self.nbolos**0.5,
        }

        self.tscache = {}

        self.tsplot = 'default'
        self.set_tsplot(**kwargs)

        print "Completed the rest of initialization in an additional %f seconds" % (time.time()-t1)

    def lookup(self, tsname):
        """
        Cache and return data...
        """
        if tsname not in self.tscache:
            t0 = time.time()
            s0 = heapy.heap().size
            print "Loading and caching %s" % tsname
            self.tscache[tsname] = self.tsplot_dict[tsname]()
            print "Loading and caching %s took %0.2g seconds and ate up %0.2g GB" % (tsname,time.time()-t0, (heapy.heap().size-s0)/1024.**3)

        return self.tscache.get(tsname)
    
    def set_tsplot(self,tsplot=None):
        """
        Options: set tsplot equal to one of these strings
        default = skysub (atmo_one-atmosphere+astrosignal)
        default_noscale (ac_bolos-atmo_one-atmosphere)
        residual (atmo_one-atmosphere-noise)
        last_astrosignal (atmo_one-atmosphere-noise+astrosignal)
        astrosignal
        dcbolos
        acbolos 
        acbolos_noscale 
        atmosphere
        default_noscale
        scale
        raw
        rawscaled
        noise
        zeromedian
        """
        if tsplot is not None:
            self.tsplot=tsplot
        if self.tsplot_dict.has_key(self.tsplot):
            self.data = self.lookup(self.tsplot) #self.tsplot_dict[self.tsplot]()
        else:
            print "No option for %s" % self.tsplot
            return
        print "Set tsplot to %s" % self.tsplot
        self._refresh()

    def readncfile(self):
          self.ncfile = netcdf.netcdf_file(self.ncfilename,'r') # NetCDF.NetCDFFile(self.ncfilename,'r')
          self.ncflags = asarray(self.ncfile.variables['flags'].data)
          self.ncbolo_params = asarray(self.ncfile.variables['bolo_params'].data)
          self.ncscans = asarray(self.ncfile.variables['scans_info'].data)
          self.timelen = self.ncflags.shape[0]
          self.scanlen = self.ncscans[0,1]-self.ncscans[0,0]
          self.whscan = asarray([arange(self.scanlen)+i for i in self.ncscans[:,0]]).ravel()
          self.nbolos = self.ncflags.shape[1]
          self.bolo_indices = asarray(nonzero(self.ncbolo_params[:,0].ravel())).ravel()
          self.nscans = self.ncscans.shape[0]
          ft = self.ncflags[self.whscan,:]
          self.ngoodbolos = self.bolo_indices.shape[0]
          self.flags = reshape(ft[:,self.bolo_indices],[self.nscans,self.scanlen,self.ngoodbolos])

    def make_noisemaps(self,save=False):
        """
        Test a variety of noisemap computations
        """
        t0=time.time()
        self.residualmap = self.mapstr['RESIDMAP'][0]
        self.weightmap = self.mapstr['WT_MAP'][0]
        self.nhitsmap = self.mapstr['NHITSMAP'][0]
        self.residsquaremap = drizzle(self.tstomap,self.noise**2,self.map.shape,self.weight*(True-self.flags))
        self.weightsquaremap = drizzle(self.tstomap,self.weight**2,self.map.shape,1.0)
        self.varscalemap = self.weightmap / (self.weightmap**2 - self.weightsquaremap)
        self.varscalemap[abs(self.weightmap**2 - self.weightsquaremap) < self.weightsquaremap/1e6] = 0
        self.rmssamplemean = (self.varscalemap * self.residsquaremap)**0.5
        self.rootresidsquaremap = self.residsquaremap**0.5
        self.smoothresid = smooth(self.residualmap,10/2.35,ignore_nan=True)
        self.smoothnoisemap = numpy.sqrt( smooth((self.residualmap-self.smoothresid)**2,10/2.35,ignore_nan=True) )
        for arrname in ("residualmap", "weightmap", "nhitsmap", "residsquaremap", "weightsquaremap", "varscalemap", "rmssamplemean", "rootresidsquaremap",):
            self.__dict__[arrname] = nantomask(self.__dict__[arrname])
            print "%20s mu=%8.2g std=%8.2g" % (arrname,self.__dict__[arrname].mean(),self.__dict__[arrname].std())
        print "Took %0.1f seconds to compute noisemaps" % (time.time()-t0)

    def save_noisemaps(self,clobber=True):
        prefix = self.filename.replace("_postiter.sav","")
        F = pyfits.open(prefix+"_map00.fits")
        for arrname in ( "residsquaremap", "weightsquaremap", "varscalemap", "rmssamplemean", "rootresidsquaremap",):
            F[0].data = self.__dict__[arrname]
            F.writeto(prefix+"_"+arrname+".fits",clobber=clobber)

    def showmap(self,colormap=cm.spectral,vmin=None,vmax=None,fignum=0,axlims=None):
      self.mapfig=figure(fignum); clf(); 
      OK = self.map == self.map
      if vmax is None:
          vmax = self.map[OK].mean()+7*self.map[OK].std()
      elif vmax=='max':
          vmax = self.map[OK].max()
      if vmin is None:
          vmin = self.map[OK].mean()-2*self.map[OK].std()
      elif vmin=='min':
          vmin = self.map[OK].min()
      self.mapim = pylab.imshow(self.map,
              vmin=vmin,vmax=vmax,
              interpolation='nearest',
              cmap=colormap); 
      self.mapim.axes.patch.set_fc('gray')
      if axlims:
          self.mapim.axes.axis(axlims)
      colorbar()
      try:
          disconnect(self.MtoT)
          disconnect(self.MtoTkey)
      except:
          pass
      self.MtoT = connect('button_press_event',self.mapclick)
      self.MtoTkey = connect('key_press_event',self.mapkeypress)
      self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)
      self.mapconnections.append(self.MtoT)
      self.mapconnections.append(self.MtoTkey)

    def showmodel(self,colormap=cm.spectral,vmin=None,vmax=None,fignum=6):
      self.modelfig=figure(fignum); clf(); 
      if vmax is None:
          vmax = self.model.mean()+7*self.model.std()
      elif vmax=='max':
          vmax = self.model.max()
      if vmin is None:
          vmin = self.model.mean()-2*self.model.std()
      elif vmin=='min':
          vmin = self.model.min()
      self.modelim = pylab.imshow(self.model,
              vmin=vmin,vmax=vmax,
              interpolation='nearest',
              cmap=colormap); 
      self.mapim.axes.patch.set_fc('gray')
      colorbar()
      try:
          disconnect(self.MtoT)
          disconnect(self.MtoTkey)
      except:
          pass
      self.MtoT = connect('button_press_event',self.mapclick)
      self.MtoTkey = connect('key_press_event',self.mapkeypress)
      self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)
      self.mapconnections.append(self.MtoT)
      self.mapconnections.append(self.MtoTkey)

    def footprint(self,tsx,tsy,scatter=False):
      mappoints = asarray(self.tstomap[self.scannum,tsy,:])

      x,y = mappoints / self.map.shape[1],mappoints % self.map.shape[1]

      if scatter:
          self.plotfig=figure(4)
          self.plotfig.clf()
          self.plotmanager = pylab.get_current_fig_manager()
          downsample_factor = 2.
          try:
              vals = self.data.data[self.scannum,tsy,:].ravel()
          except TypeError:
              vals = self.data[self.scannum,tsy,:].ravel()
          try:
              flags = self.flags.data[self.scannum,tsy,:].ravel()
          except TypeError:
              flags = self.flags[self.scannum,tsy,:].ravel()
          flagvals = vals*(flags==0)
          self.footim = pylab.imshow(gridmap(x,y,flagvals,downsample_factor=downsample_factor,xsize=72,ysize=72)
                  ,interpolation='bilinear',extent=[0,36*7.2,0,36*7.2])
          pylab.xlabel('Arcseconds')
          pylab.ylabel('Arcseconds')
          self.footcb = pylab.colorbar()
          try:
              self.footscatter = pylab.scatter(7.2*(x-min(x))/downsample_factor,7.2*(y-min(y))/downsample_factor,c=self.data[self.scannum,tsy,:],s=40)

          except TypeError:
              self.footscatter = pylab.scatter(7.2*(x-min(x))/downsample_factor,7.2*(y-min(y))/downsample_factor,c=self.data.data[self.scannum,tsy,:],s=40)

      else:
          try:
              self.fp2[0].set_visible(False)
              self.fp3[0].set_visible(False)
              self._refresh()
          except:
              pass

          figure(0)
          myaxis = self.mapfig.axes[0].axis()
          self.fp3 = plot(y,x,'ro')
      #    self.fp1 = plot(y,x,'b+')
          self.fp2 = plot(y,x,'wx')
          self.mapfig.axes[0].axis(myaxis)
      self._refresh()

    def bolomap(self,bolonum):
        bolonum = numpy.round(bolonum)
        self.bolofig = pylab.figure(5)
        self.bolofig.clf()
        self.bolommap = numpy.zeros(self.map.shape)
        self.bolonhits = numpy.zeros(self.map.shape)
        self.bolommap.flat[self.tstomap[:,:,bolonum].ravel()] += self.data[:,:,bolonum].ravel()
        # only add hits for non-zero, non-NAN elements
        hits = (self.data[:,:,bolonum].ravel() != 0) * (self.data[:,:,bolonum].ravel() == self.data[:,:,bolonum].ravel())
        self.bolonhits.flat[self.tstomap[:,:,bolonum].ravel()] += hits
        self.bolommap /= self.bolonhits
        self.bolomapim = pylab.imshow(self.bolommap,interpolation='nearest',origin='lower')
        title("Bolometer %i" % bolonum)
        try:
            disconnect(self.gfit_connection)
        except:
            pass
        self.gfit_connection = connect('button_press_event',self.gfit_bolomap)
        pylab.colorbar()
        return self.bolommap

    def gfit_map(self,event,map,ax=None):
        if ax is None:
            if self.bolofig is not None:
                tb = self.bolofig.canvas.manager.toolbar
            elif self.mapfig is not None:
                tb = self.mapfig.canvas.manager.toolbar
            else:
                raise Exception("Can't gaussian fit - no canvas exists?")
        else: 
            tb = ax.figure.canvas.manager.toolbar
        if tb.mode=='' and event.xdata and event.ydata:
            errmask = 1e10 * getmask(map)+1.0
            gf = gaussfitter.gaussfit(masktozero(map),
                    err=errmask,
                    params=   [0,0,event.xdata,event.ydata,0,0,0],
                    usemoment=[1,1,          0,          0,1,1,1],
                    rotate=1)
            self.gaussfit = gf
            if ax:
                ax.add_patch(Ellipse(gf[2:4],numpy.max(gf[4:6])*2.35,numpy.min(gf[4:6])*2.35,gf[6],fill=False))
            print "Guess: %i,%i" % (event.xdata,event.ydata)," Fit peak: %g  Background: %g  X,Y position: %f,%f  X,Y FWHM: %f,%f   Angle: %f" % (gf[1],gf[0],gf[2],gf[3],gf[4]*2.35*7.2,gf[5]*2.35*7.2,gf[6])

    def gfit_bolomap(self,event):
        self.gfit_map(event,self.bolommap,self.bolomapim.axes)
        self.bolofig.canvas.draw()

    def fiteachbolo(self):
        self.gfits = numpy.array([ gaussfitter.gaussfit(masktozero(self.bolomap(bolonum)),return_all=1) for bolonum in range(self.nbolos)])
        self.mapfit = gaussfitter.gaussfit(self.map,return_all=1)
        self.distmapfig = pylab.figure(7)
        self.distmapplot= pylab.errorbar(self.gfits[:,0,2]-self.mapfit[0,2],
                self.gfits[:,0,3]-self.mapfit[0,3],
                self.gfits[:,1,2],self.gfits[:,1,3],fmt=None,capsize=3)
        return self.gfits

    def footmovie(self,y1,y2,movie=False,moviedir='scanmovie/',logscale=False,dosmooth=True):
        y1 = numpy.round(y1)
        y2 = numpy.round(y2)
        print "Making footprint movie from %i to %i" % (y1,y2)
        self.footscatter.set_visible(False)
        #if isinstance(self.data,numpy.ma.masked_array):
        #    plotdata = self.data.data
        #    allflags = self.flags.data
        #else:
        plotdata = self.data
        allflags = self.flags


        mappoints = asarray(self.tstomap[self.scannum,y1:y2,:])
        x,y = mappoints / self.map.shape[1],mappoints % self.map.shape[1]
        vals = numpy.reshape( plotdata[self.scannum,y1:y2,:], [abs(y2-y1),plotdata.shape[2]] )
        flags = numpy.reshape( allflags[self.scannum,y1:y2,:], [abs(y2-y1),plotdata.shape[2]] )
        flagvals = vals*(flags==0)
        mapcube = array([ gridmap(x[ii,:],y[ii,:],flagvals[ii,:],downsample_factor=2,xsize=72,ysize=72,dosmooth=dosmooth)
            for ii in xrange(y2-y1) ])
        mapcube[mapcube!=mapcube] = 0

        if logscale:
            #mapcube-=mapcube.min()
            mapcube=numpy.arcsinh(mapcube/logscale)
            #self.footim.set_norm(matplotlib.colors.LogNorm())
        vmin = mapcube.min()
        vmax = mapcube.max()
        #self.footim.set_clim((vmin,vmax))
        #self.footcb.set_clim((vmin,vmax))
        #self.footim.set_array(mapcube[0,:,:])
        figure(4)
        self.footim = imshow(mapcube[0,:,:],vmin=vmin,vmax=vmax)

        self.footcb.ax.clear()
        self.footcb = pylab.colorbar(cax=self.footcb.ax) #self.footim,cax=self.footcb.ax,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax))

        pylab.draw()
        for tsy in xrange(y2-y1):
            self.footim.set_array(mapcube[tsy,:,:])
            #self.plotfig.show()
            #pylab.draw()
            self.plotmanager.canvas.draw()
            if movie:
                pylab.savefig(moviedir+'%04i.png' % tsy)
        
        return mapcube

 
    def set_plotscan_data(self,scannum,data=None,flag=True):
        if data is not None and flag:
            self.plane = data *(self.flags[scannum,:,:] == 0)
        elif data is not None:
            self.plane = data
        elif flag:
            self.plane = self.data[scannum,:,:] * (self.flags[scannum,:,:]==0)
        else:
            self.plane = self.data[scannum,:,:] 

    def plotscan(self, scannum, fignum=1, button=1, data=None, flag=True, logscale=False):
      if self.connected:
          self.dcon()
      self.scannum = scannum
      self.set_plotscan_data(scannum,flag=flag,data=data)
      self.fignum = fignum
      if logscale:
          plotdata = log10(abs(self.plane)) * sign(self.plane)
      else:
          plotdata = self.plane
      if self.scanim is None:
          self.flagfig = figure(fignum+1,figsize=[16,12]); clf()
          if texOn:
            self.flagtitle = pylab.title("Flags for Scan "+str(self.scannum)+" in "+self.ncfilename.replace("_","\\_"));
          else:
            self.flagtitle = pylab.title("Flags for Scan "+str(self.scannum)+" in "+self.ncfilename);
          xlabel('Bolometer number'); ylabel('Time (0.%0.2fs)' % self.sample_interval)
          self.flagim = pylab.imshow(self.flags[scannum,:,:],interpolation='nearest',
                  origin='lower',aspect=self.aspect)
          self.flagcb = pylab.colorbar()
          self.datafig = figure(fignum,figsize=[16,12]);clf();
          if texOn:
            self.datatitle = pylab.title("Scan "+str(self.scannum)+" in "+self.ncfilename.replace("_","\\_"));
          else:
            self.datatitle = pylab.title("Scan "+str(self.scannum)+" in "+self.ncfilename);
          xlabel('Bolometer number'); ylabel('Time (0.%0.2fs)' % self.sample_interval)
          self.dataim = pylab.imshow(plotdata,interpolation='nearest',
                  origin='lower',aspect=self.aspect)
          self.datacb = pylab.colorbar()
      else:
          self.flagim.set_array(self.flags[scannum,:,:])
          self.flagcb = pylab.colorbar(cax=self.flagcb.ax)
          self.dataim.set_array(plotdata)
          self.datacb = pylab.colorbar(cax=self.datacb.ax)
          if texOn:
            self.flagtitle.set_text("Flags for Scan "+str(self.scannum)+" in "+self.ncfilename.replace("_","\\_"))
            self.datatitle.set_text("Scan "+str(self.scannum)+" in "+self.ncfilename.replace("_","\\_"))
          else:
            self.flagtitle.set_text("Flags for Scan "+str(self.scannum)+" in "+self.ncfilename)
            self.datatitle.set_text("Scan "+str(self.scannum)+" in "+self.ncfilename)
      self.showrects()
      self.showlines()
      self.cursor = Cursor(self.dataim.axes,useblit=True,color='black',linewidth=1)
      self._refresh()
      self.reconnect()
      """
      self.md  = connect('button_press_event',self.mouse_down_event)
      self.mu  = connect('button_release_event',self.mouse_up_event)
      self.key = connect('key_press_event',self.keypress)
      self.connections.append(self.md)
      self.connections.append(self.mu)
      self.connections.append(self.key)
      """

    def flagpoint(self, i, j, button):
        if button==1:
            self.flags[self.scannum,round(j),round(i)] += 1
#          print >>self.outfile,\
#                  "flag_manual,'%s',bad_bolos=[%i],bad_time=[[%i,%i]],/doboth"\
#                  % (self.ncfilename,i,self.scannum,j)
        elif button==2:
            self.flags[self.scannum,round(j),round(i)] -=\
                    (self.flags[self.scannum,round(i),round(j)] > 0)
#          print >>self.outfile,\
#                  "undo_flag,'%s',bad_bolos=[%i],bad_time=[[%i,%i]],/doboth" \
#                  % (self.ncfilename,i,self.scannum,j)
        elif button==3 or button=='d':
            for p in self.rectangles[self.scannum]:
                if p.get_window_extent().contains(self.event.x,self.event.y):
                    p.set_visible(False)
                    self.rectangles[self.scannum].remove(p)
                    x1,x2 = (p.get_x()+.5*sign(p.get_width()) ,
                            p.get_x()-.5*sign(p.get_width())+p.get_width() )
                    y1,y2 = (p.get_y()+.5*sign(p.get_height()),
                            p.get_y()-.5*sign(p.get_height())+p.get_height() )
                    if y1 > y2:
                      y1,y2=y2,y1
                    if x1 > x2:
                      x1,x2=x2,x1
                    if p.get_fc() == 'black':
                        self.flags[self.scannum,y1:y2+1,x1:x2+1] -= (
                                (self.flags[self.scannum,y1:y2+1,x1:x2+1] != 0) 
                                * sign(self.flags[self.scannum,y1:y2+1,x1:x2+1]) )
                          # only subtract 1 so that overlaps aren't completely unflagged
                    elif p.get_fc() == 'blue':
#                      print x1,x2,y1,y2,p.xy,p.get_width(),p.get_height()
                        self.flags[self.scannum,y1:y2+1,x1:x2+1] *= (
                                -1*(self.flags[self.scannum,y1:y2+1,x1:x2+1] < 0) 
                                + 1*(self.flags[self.scannum,y1:y2+1,x1:x2+1] > 0) )
#                  print >>self.outfile,"Removed object with center %f,%f"\
#                          % (p.get_x(),p.get_y())
                    break
            self._refresh()
 
    def flag_box(self,x1,y1,x2,y2,button):
#      x = (x1+x2)/2.
#      y = (y1+y2)/2.
        x1i = int(round(x1))
        x2i = int(round(x2))
        y1i = int(round(y1))
        y2i = int(round(y2))
        w = (x1i-x2i)+sign(x1i-x2i)
        if abs(w) == 0:
            w = sign(x1-x2)
        h = (y1i-y2i)+sign(y1i-y2i)
        if abs(h) == 0:
            h = sign(y1-y2)
        if y1==y2:
            h = 1
        if x1==x2:
            w = 1
        x2 = x2i-.5*sign(w)
        y2 = y2i-.5*sign(h)
        x1 = x2i+w
        y1 = y2i+h
        yrange = [min(y1i,y2i),min(y1i,y2i)+abs(h)]
        xrange = [min(x1i,x2i),min(x1i,x2i)+abs(w)]
        scannum = self.scannum
        if button==1:
            self.flags[scannum,yrange[0]:yrange[1],xrange[0]:xrange[1]] += 1
            p = matplotlib.patches.Rectangle(xy=(x2,y2), width=w, height=h,
                    facecolor='black',transform=gca().transData)
            gca().add_patch(p)
            p.set_visible(True)
            p.set_alpha(.5)
#          self.axis.draw()
            self.rectangles[self.scannum].append(p)
            self._refresh()
#          print x,y,w,h,p
#          print >>self.outfile,\
#                  "flag_manual,'%s',bolorange=[%i,%i],timerange=[%i,%i],scanrange=%i" \
#                  % (self.ncfilename,x1,x2,y1,y2,self.scannum)
        elif button==2:
           # this won't work right; I need some way to make it undo-able  
           # <--- I don't know if that's true any more (11/10/08)
            unflagreg = self.flags[scannum,yrange[0]:yrange[1],xrange[0]:xrange[1]] 
            unflagreg[unflagreg > 0] = 0
            p = matplotlib.patches.Rectangle(xy=(x2,y2), width=w, height=h,
                    facecolor='blue',transform=gca().transData)
            gca().add_patch(p)
            p.set_visible(True)
            p.set_alpha(.5)
            self.rectangles[self.scannum].append(p)
            self._refresh()
#          print >>self.outfile,\
#                  "undo_flag,'%s',bolorange=[%i,%i],timerange=[%i,%i],scanrange=%i" \
#                  % (self.ncfilename,x1,x2,y1,y2,self.scannum)
        elif button=='d':
            for p in self.rectangles[self.scannum]:
                if p.get_window_extent().contains(self.event.x,self.event.y):
                    p.set_visible(False)
                    self.rectangles[self.scannum].remove(p)
                    x1,x2 = (p.get_x()+.5*sign(p.get_width()) ,p.get_x()
                            -.5*sign(p.get_width())+p.get_width() ) 
                    y1,y2 = (p.get_y()+.5*sign(p.get_height()),p.get_y()
                            -.5*sign(p.get_height())+p.get_height() )
                    if y1 > y2:
                      y1,y2=y2,y1
                    if x1 > x2:
                      x1,x2=x2,x1
                    if p.get_fc() == 'black':
                        self.flags[self.scannum,y1:y2+1,x1:x2+1] -= (
                                (self.flags[self.scannum,y1:y2+1,x1:x2+1] != 0) 
                                * sign(self.flags[self.scannum,y1:y2+1,x1:x2+1]) )
                          # only subtract 1 so that overlaps aren't completely unflagged
                    elif p.get_fc() == 'blue':
#                      print x1,x2,y1,y2,p.xy,p.get_width(),p.get_height()
                        self.flags[self.scannum,y1:y2+1,x1:x2+1] *= (
                                -1*(self.flags[self.scannum,y1:y2+1,x1:x2+1] < 0) 
                                + 1*(self.flags[self.scannum,y1:y2+1,x1:x2+1] > 0) )
#                  print >>self.outfile,"Removed object with center %f,%f"\
#                          % (p.get_x(),p.get_y())
                    break
            self._refresh()
        elif button==3:
            self.maparrow(x2i,y2i)
 
    def flag_bolo(self,x,button):
        if button=='b' and self.PCAflag:
            x=round(x)
            PCAsub = self.data[self.scannum,:,:] - outer(self.plane[:,x],ones(self.plane.shape[1]))
            self.PCAflag = False
            self.plotscan(self.scannum,data=PCAsub)
        elif button=='b' and not self.PCAflag:
#          print x,round(x),button,self.data.shape
            x=round(x)
            h=self.data.shape[1]
            self.flags[self.scannum,0:h,x] += 1
            p = matplotlib.lines.Line2D([x,x],[0,h],\
                    color='black',transform=gca().transData)
            gca().add_line(p)
            p.set_visible(True)
            self.lines[self.scannum].append(p)
            self._refresh()
#          print >>self.outfile,\
#                  "flag_manual,'%s',bolorange=[%i],scanrange=%i" \
#                  % (self.ncfilename,x,self.scannum)

    def flag_times(self,y1,y2):
        y1=max([round(y1),0])
        y2=min([round(y2),self.data.shape[1]])
        w=self.data.shape[2]
        self.flags[self.scannum,y1:y2,0:w] += 1
        p1 = matplotlib.lines.Line2D([0,w],[y1,y1],\
                color='black',transform=gca().transData)
        p2 = matplotlib.lines.Line2D([0,w],[y2,y2],\
                color='black',transform=gca().transData)
        gca().add_line(p1)
        gca().add_line(p2)
        p1.set_visible(True)
        p2.set_visible(True)
        self.lines[self.scannum].append(p1)
        self.lines[self.scannum].append(p2)
        self._refresh()

    def unflag_times(self,y,button):
        y1=max([round(y1),0])
        y2=min([round(y2),self.data.shape[1]])
        w=self.data.shape[2]
        flagarea = self.flags[self.scannum,y1:y2,0:w]
        flagarea[flagarea>0] -= 1
        for y in (y1,y2):
            for l in self.lines[self.scannum]:
                if l._y[0] == y and l._y[1] == y:
                    l.set_visible(False)
                    self.lines[self.scannum].remove(l)
        self._refresh()

    def unflag_bolo(self,x,button):
        if button=='B':
#          print x,round(x),button,self.data.shape
            x=round(x)
            h=self.data.shape[1]
            for l in self.lines[self.scannum]:
                if l._x[0] == x and l._x[1] == x:
                    self.flags[self.scannum,0:h,x] -= 1
                    l.set_visible(False)
                    self.lines[self.scannum].remove(l)
                    self._refresh()
            if self.flags[self.scannum,0:h,x].max() > 0:
                arr = self.flags[self.scannum,0:h,x]
                arr[arr>0] -= 1

    def plot_column(self,tsx,clear=True,timestream='data', color='k', fignum=7, **kwargs):
        self.bolofig=figure(fignum)
        if clear: self.bolofig.clf()
        if self.powerspec_plotted:
            xlen = self.plane.shape[0]
            pylab.plot(fftfreq(xlen,d=self.sample_interval)[:xlen/2],self.plane[:xlen/2,numpy.round(tsx)],
                    linewidth=0.5,color=color)
            xlabel("Frequency (Hz)")
            ylabel("Power (Jy$^2$)")
            title("Bolo %i" % round(tsx))
        else:
            title("Bolo %i" % round(tsx))
            xlabel("Time Samples")
            ylabel("Flux (Jy)")
            if timestream == 'data':
                pylab.plot(self.plane[:,numpy.round(tsx)],
                        linewidth=0.5,color=color,**kwargs)
            elif timestream in self.tsplot_dict.keys():
                plane = self.lookup(timestream)[self.scannum,:,:]
                pylab.plot(plane[:,numpy.round(tsx)],
                        linewidth=0.5,color=color,**kwargs)

    def plot_line(self,tsy,clear=True,fignum=7):
        self.bolofig=figure(fignum)
        if clear: self.bolofig.clf()
        if self.PCAflag:
            pylab.plot(self.efuncarr[numpy.round(tsy),:])
        else:
            title("Line %i" % round(tsy))
            xlabel("Bolometer")
            ylabel("Flux (Jy)")
            pylab.plot(self.plane[numpy.round(tsy),:])

    def dcon(self):
        self.connected = False
        disconnect(self.md)
        disconnect(self.mu)
        disconnect(self.key)
        disconnect(self.MtoT)
        for i in self.connections:
            self.mapfig.canvas.mpl_disconnect(i)
            try:
                self.datafig.canvas.mpl_disconnect(i)
            except:
                continue
            try:
                self.flagfig.canvas.mpl_disconnect(i)
            except:
                continue
        for i in self.mapconnections:
            self.mapfig.canvas.mpl_disconnect(i)

    def reconnect(self):
      self.connected = True
      self.MtoT = self.mapfig.canvas.mpl_connect('button_press_event',self.mapclick)
      self.MtoTkey = self.mapfig.canvas.mpl_connect('key_press_event',self.mapkeypress)
      self.md  = self.datafig.canvas.mpl_connect('button_press_event',self.mouse_down_event)
      self.mu  = self.datafig.canvas.mpl_connect('button_release_event',self.mouse_up_event)
      self.key = self.datafig.canvas.mpl_connect('key_press_event',self.keypress)
      self.connections.append(self.md)
      self.connections.append(self.mu)
      self.connections.append(self.key)
      self.mapconnections.append(self.MtoT)
      self.mapconnections.append(self.MtoTkey)

    def close(self,write=True):
        """ close the ncdf file and the graphics windows
        and flush everything to file"""
        if self.open == 1:
             self.open = 0
             self.dcon()
             if write and self.ncfile:
                 self.write_ncdf()
             elif write:
                self.writeflags()
#           self.outfile.close()
             if self.mapfig is not None:
                 self.mapfig.clf()
             if self.datafig is not None:
                 self.datafig.clf()
             if self.bolofig is not None:
                 self.bolofig.clf()
             if self.plotfig is not None:
                 self.plotfig.clf()
             if self.flagfig is not None:
                 self.flagfig.clf()

    def writeflags(self):
        flagged_scans = np.array([ii for ii in xrange(self.nscans) if all(self.flags[ii,:,:])])
        flagged_bolos = np.array([ii for ii in xrange(self.nbolos) if all(self.flags[:,:,ii])])
        flagged_scans.tofile(self.pathprefix+self.fileprefix+"_flagged_scans.txt"," ")
        flagged_bolos.tofile(self.pathprefix+self.fileprefix+"_flagged_bolos.txt"," ")
        tempdata = self.data
        if self.tsfile:
            self.tsfile[0].data = asarray(self.flags,dtype='int')
            self.tsfile.writeto(self.flagfn,clobber=True)
            self.tsfile[0].data = tempdata
        elif self.header:
            self.flagfits = pyfits.PrimaryHDU(asarray(self.flags,dtype='int'),self.header)
            self.flagfits.writeto(self.flagfn,clobber=True)
            print "flags_to_ncdf,'%s','%s'" % (self.flagfn,self.filename)
            self.idl_writeflags()

    def idl_writeflags(self):
        flagcmd = "flags_to_ncdf,'%s','%s'" % (self.flagfn,self.filename)
        idlcmd = "/Applications/itt/idl/idl/bin/idl"
        cmd = '%s -e "%s"' % (idlcmd,flagcmd)
        os.environ.update({"IDL_STARTUP":"/Users/adam/work/bolocam/.idl_startup_bgps.pro"})
        print cmd
        output = os.popen(cmd)
        print "".join(output.readlines())
        output.close()

    def mapclick(self,event):
        if event.xdata == None:
            return
        clickX = round(event.xdata)
        clickY = round(event.ydata)
        if event.button == 1:
            self.tsarrow(clickX,clickY)
        elif event.button == 2:
            self.find_all_points(clickX,clickY)
        elif event.button == 3:
            self.hist_all_points(clickX,clickY)

    def mapkeypress(self,event):
        if event.inaxes is None: return
        else:
            clickX = round(event.xdata)
            clickY = round(event.ydata)
        if event.key == 'c':
            self.toggle_currentscan()
        elif event.key == 'G':
            self.gfit_map(event,self.map,self.mapim.axes)
            self.mapfig.canvas.draw()
        elif event.key == '.':
            if event.xdata == None:
                return
            self.find_all_points(clickX,clickY)
            self.tsarrow(clickX,clickY)
        elif event.key == "r":
            self.showmap()
        elif event.key == 'a':
            self.tsarrow(clickX,clickY)
        elif event.key == 'm':
            self.find_all_points(clickX,clickY)
        elif event.key == 'h':
            self.hist_all_points(clickX,clickY)

    def keypress(self,event):
        set_lastkey=True
        if event.inaxes is None: return
        elif event.key == 'n':
            if self.scannum < self.maxscan-1:
                self.plotscan(self.scannum+1)
            else:
                print "At last scan, can't go further"
        elif event.key == 'p' or event.key == 'N':
            if self.scannum > 0:
                self.plotscan(self.scannum-1)
            else:
                print "At first scan, can't go further back"
        elif event.key == 'e': 
            self.expsub()
        elif event.key == 'P': # PCA
            self.plotscan(self.scannum,data=efuncs(self.plane),flag=False,logscale=True)
            self.PCAflag = True
        elif event.key == 'q':
            self.close()
        elif event.key == 'Q':
            self.close(write=False)
        elif event.key == '.':
            self.maparrow(round(event.xdata),round(event.ydata))
        elif event.key == 'f':
            self.footprint(round(event.xdata),round(event.ydata))
        elif event.key == 'F':
            self.footprint(round(event.xdata),round(event.ydata),scatter=True)
        elif event.key == 'R': # reverse order of boxes
            self.rectangles[self.scannum].reverse()
        elif event.key == 'r': # redraw
            self.plotscan(self.scannum)
        elif event.key == 'M': # flag highest point
            self.flags[self.scannum,:,:].flat[self.plane.argmax()] += 1
            self.plane.flat[self.plane.argmax()] = 0
        elif event.key == 'm': # flag lowest point
            self.flags[self.scannum,:,:].flat[self.plane.argmin()] += 1
            self.plane.flat[self.plane.argmin()] = 0
        elif event.key == 'd':
            self.flag_box(self.x1,self.y1,self.x2,self.y2,'d')
        elif event.key == 't' or event.key == 'T':
            if self._lastkey == 't' or self._lastkey == 'T':
                self._y2 = numpy.ceil(event.ydata)
                if event.key == 'T':
                    self.unflag_times(self._y1,self._y2)
                elif event.key =='t':
                    self.flag_times(self._y1,self._y2)
                self._lastkey = None
                set_lastkey = False
            else:
                self._y1 = numpy.floor(event.ydata)
        elif event.key == 's' or event.key == 'w': # "whole" scan
            self.flags[self.scannum,:,:] += 1
        elif event.key == 'S' or event.key == 'W':
            self.flags[self.scannum,:,:] -= (self.flags[self.scannum,:,:] > 0)
        elif event.key == 'b':
            self.flag_bolo(event.xdata,event.key)
        elif event.key == 'B':
            self.unflag_bolo(event.xdata,event.key)
        elif event.key == 'c':
            self.toggle_currentscan()
        elif event.key == 'C':
            self.plot_column(event.xdata)
        elif event.key == 'L':
            self.plot_line(event.ydata)
        elif event.key == 'z':
            self.powerspec()
        elif event.key == 'Z':
            self.powerspec_whole(event.xdata)
        elif event.key == 'j':
            self.timestream_whole(event.xdata)
        elif event.key == 'a':
            if self._lastkey == 'a':
                self._y2 = round(event.ydata)
                self.skymovie = self.footmovie(self._y1,self._y2)
                self._lastkey = None
                set_lastkey = False
            else:
                self._y1 = round(event.ydata)
                self.footprint(round(event.xdata),round(event.ydata),scatter=True)
        elif event.key == 'o':
            self.bolomap(event.xdata)
        elif event.key == 'v':
            x,y = round(event.xdata),round(event.ydata)
            vpt = self.data[self.scannum,y,x]
            fpt = self.flags[self.scannum,y,x]
            xmap = self.tstomap[self.scannum,y,x] / self.map.shape[1]
            ymap = self.tstomap[self.scannum,y,x] % self.map.shape[1]
            print "Value at %i,%i: %f  Flagged=%i  Maps to: %i,%i" % (x,y,vpt,fpt,xmap,ymap)
        elif event.key == '?':
            print self.help
        if set_lastkey: 
            self._lastkey = event.key

    def find_all_points(self,x,y):
        mappoint = y * self.map.shape[1] + x
        self.timepoints =  nonzero(self.tstomap == mappoint)
        wtavg = (self.mapped_timestream[self.timepoints]*self.weight[self.timepoints]).sum() / self.weight[self.timepoints].sum()
        # not a real thing wtsclavg = (self.mapped_timestream[self.timepoints]*self.weight[self.timepoints]*self.scalearr[self.timepoints]).sum() / (self.weight[self.timepoints]*self.scalearr[self.timepoints]).sum()
        uwtavg = self.mapped_timestream[self.timepoints].mean()
        medavg = median(self.mapped_timestream[self.timepoints])
        Hmad = MAD(self.mapped_timestream[self.timepoints])
        Hstd = std(self.mapped_timestream[self.timepoints])
        print ""
        print "Location: %i,%i" % (x,y)
        print "Map value: %f   Weighted average: %f   Unweighted Average: %f  Median: %f" % (self.map[y,x],wtavg,uwtavg,medavg)
        print "MAD: %f   StdDev: %f" % (Hmad,Hstd)
        print "scan,bolo,time: %12s%12s%12s%12s%12s%12s%12s%12s" % ('mapped','mapped_astr','astro','noise','residual','flags','weight','scale')
        for ii,jj,kk in transpose(self.timepoints):
            if self.flags[ii,jj,kk]:
                print "%4i,%4i,%4i: %12s%12s%12s%12f%12s%12i%12s%12s" % (ii,kk,jj,"","","",self.noisemap[y,x],"",self.flags[ii,jj,kk],"","")
            else:
                print "%4i,%4i,%4i: %12f%12f%12f%12f%12f%12i%12f%12f" % (ii,kk,jj,
                    self.mapped_timestream[ii,jj,kk],
                    self.mapped_astrosignal[ii,jj,kk],
                    self.astrosignal[ii,jj,kk],
                    self.noisemap[y,x],
                    self.noise[ii,jj,kk],
                    self.flags[ii,jj,kk],
                    self.weight[ii,jj,kk],
                    self.scalearr[ii,jj,kk])

    def expsub(self):
        for ii in xrange(self.nbolos):
            if hasattr(self.plane,'mask'):
                if self.plane.mask[:,ii].sum() < self.plane.shape[0] - 2:
                    self.plane[:,ii] = expsub_line(self.plane[:,ii])
        self.plotscan(self.scannum, data=self.plane, flag=False)

    def hist_all_points(self,x,y,clear=True,timestream='mapped_timestream'):
        mappoint = y * self.map.shape[1] + x
        self.timepoints =  nonzero(self.tstomap == mappoint)
        if timestream in self.tsplot_dict.keys():
            TS = self.lookup(timestream) #tsplot_dict[timestream]()
        else:
            raise KeyError("Timestream %s is not valid" % (timestream))
        wtavg = (TS[self.timepoints]*self.weight[self.timepoints]).sum() / self.weight[self.timepoints].sum()
        uwtavg = TS[self.timepoints].mean()
        medavg = median(TS[self.timepoints])
        Hmad = MAD(TS[self.timepoints])
        Hstd = std(TS[self.timepoints])
        datapts = TS[self.tstomap==mappoint]
        self.plotfig=figure(4)
        self.plotfig.clear()
        OK  = asarray(datapts == datapts)
        if hasattr(datapts,'mask'):
            OK *= (datapts.mask == False)
        n,bins,patches = hist(asarray(datapts[OK]),histtype='step',color='k',linewidth=2)
        vlines(wtavg,0,max(n),color='k',linestyles=':',label="Weighted: %0.4g $\\pm$ %0.4g" % (wtavg,Hstd))
        #vlines(wtavg,0,max(n),color='k',linestyles=':',label="Std: %0.4g" % Hstd)
        fill_betweenx([0,max(n)],[wtavg-Hstd]*2,[wtavg+Hstd]*2,color='k',alpha=0.1,label="Std: %0.4g" % Hstd)
        vlines(uwtavg,0,max(n),color='b',linestyles='-.',label="Unweighted: %0.4g" % uwtavg)
        vlines(medavg,0,max(n),color='g',linestyles='--',label="Median: %0.4g $\\pm$ %0.4g" % (medavg,Hmad))
        fill_betweenx([0,max(n)],[medavg-Hmad]*2,[medavg+Hmad]*2,color='g',alpha=0.1,label="MAD: %0.4g" % Hmad)
        vlines(self.model[y,x],0,max(n),color='purple',linestyles='--',label="Model: %0.4g" % Hmad)
        Ctemp = matplotlib.collections.CircleCollection([0],facecolors='k',edgecolors='k')
        Ctemp.set_label('Map Value: %0.4g' % (self.map[y,x]))
        self.plotfig.axes[0].add_collection(Ctemp)
        L=legend(loc='best')
        L.draggable(True)
        title("%s pixel %i,%i" % (self.filename,x,y))
        xlabel('Flux (Jy or Volts)')

    def tsarrow(self,x,y):
        if self.debug: print "tsarrow at %f,%f" % (x,y)
        #      xy = [clickX,clickY]

        # this took a little thinking:
        # the Y axis has HUGE variation, X has small....
        mappoint = y * self.map.shape[1] + x
        self.timepoints =  nonzero(self.tstomap == mappoint)

        matchpts = list(nonzero(self.timepoints[0] == self.scannum))[0]

#      print mappoint,clickX,clickY,self.timepoints,outer(xy,self.map.shape)
#      for i in outer(xy,self.map.shape).ravel():
#          print i," :  ",nonzero(self.tstomap==mappoint)

#      print matchpts,mappoint,self.timepoints

        if self.connected:
            for a in self.arrows:
                a.set_visible(False)
            for a in self.arrows:
                self.arrows.remove(a)
            for i in list(matchpts):
                if self.debug: print "i shape: ",i.shape, " matchpts ",matchpts
                i = int(i)
                t,b = self.timepoints[1][i],self.timepoints[2][i]
#              print "T,b,i  ",t,b,i
#              print "Does t = []?",t == []
#              print "Is t >= 0?",t >= 0
#              arrow = FancyArrow(t-5,b-5,5,5)
#              self.datafig.axes[0].add_patch(arrow)
                figure(self.fignum)
                ax = self.datafig.axes[0]
                # redundant? self.datafig.sca(self.datafig.axes[0])
                #arrow = self.datafig.axes[0].arrow(t-5,b-5,5,5)
                a1 = ax.arrow(b-3,t-3,6,6,head_width=0,facecolor='black')
                a2 = ax.arrow(b-3,t+3,6,-6,head_width=0,facecolor='black')
                a1.set_visible(True)
                a2.set_visible(True)
#              print a,t,b
                self.arrows.append(a1)
                self.arrows.append(a2)
            self._refresh()

    def maparrow(self,tsx,tsy):

#    scanpoint = self.scannum*self.flags.shape[1]*self.flags.shape[2]\
          #    + y*self.flags.shape[0] + x
#    print tsx,tsy
      mappoint = self.tstomap[self.scannum,tsy,tsx]
      x,y = mappoint / self.map.shape[1],mappoint % self.map.shape[1]

      for a in self.maparrows:
          a.set_visible(False)
      for a in self.maparrows:
          self.maparrows.remove(a)
      figure(0)
      ax = self.mapfig.axes[0]
      a1 = ax.arrow(y+2,x+2,-4,-4,head_width=0,facecolor='black',
              length_includes_head=True,head_starts_at_zero=False)
      a2 = ax.arrow(y-2,x+2,4,-4,head_width=0,facecolor='black',
              length_includes_head=True,head_starts_at_zero=False)
      a1.set_visible(True)
      a2.set_visible(True)
      self.maparrows.append(a1)
      self.maparrows.append(a2)
      self._refresh()

    def toggle_currentscan(self):
        if self.currentscan == 0:
            xarr = self.tstomap[self.scannum,:,:] / self.map.shape[1]
            yarr = self.tstomap[self.scannum,:,:] % self.map.shape[1]
            x0,x1 = xarr.min(),xarr.max()
            y0,y1 = yarr.min(),yarr.max()
            self.mapfig.axes[0].axis([y0,y1,x0,x1])
            self.currentscan = 1
            self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)
        elif self.currentscan == 1:
            self.mapfig.axes[0].axis([0,self.map.shape[1],0,self.map.shape[0]])
            self.currentscan = 0
            self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)


    def showrects(self):
        ax = gca()
        for p in self.rectangles[self.scannum]:
            p.set_transform(ax.transData)
            ax.add_patch(p)

    def showlines(self):
        ax = gca()
        for l in self.lines[self.scannum]:
            l.set_transform(ax.transData)
            ax.add_line(l)

    def reset(self):
      """ Reset flags after the update function is called.
          Mouse is tracked separately.
          """
      self.limits_changed = 0
      self.got_draw = False

    def mouse_up_event(self, event):
      if event.inaxes is None: return
      self.mouse_up = True
      self.x2 = event.xdata
      self.y2 = event.ydata
      self.event = event
      tb = get_current_fig_manager().toolbar
      if tb.mode=='' and not self.PCAflag:
          self.flag_box(self.x1,self.y1,self.x2,self.y2,event.button)
#        if abs(self.x2-self.x1) > 1 or abs(self.y2-self.y1) > 1:
#        else:
#            self.flagpoint(self.x1,self.y1,event.button)

    def mouse_down_event(self, event):
      if event.inaxes is None: return
      self.mouse_up = False
      self.x1 = event.xdata
      self.y1 = event.ydata

    def powerspec(self):
        self.powerspectra = real(fft(masktozero(self.plane),axis=0) * conj(fft(masktozero(self.plane),axis=0)))
        self.plotscan(self.scannum,data=self.powerspectra,flag=False,logscale=True)
        ylabel('Frequency')
        self.powerspec_plotted = True

    def powerspec_whole(self,bolonum=0,recompute=False,timestream='data',clear=True,fignum=4,logx=False,logy=True,color='k'):
        if self.powerspectra_whole is None or recompute:
            if timestream == 'data':
                wholedata = reshape(self.data,[self.data.shape[0]*self.data.shape[1],self.data.shape[2]])
            elif self.tsplot_dict.has_key(timestream):
                data = self.lookup(timestream) #tsplot_dict[timestream]()
                wholedata = reshape(data,[data.shape[0]*data.shape[1],data.shape[2]])
            else:
                raise KeyError("Timestream %s is not valid." % timestream)
            if hasattr(wholedata,'data'): wholedata = wholedata.data
            wholedata[wholedata!=wholedata] = 0
            self.powerspectra_whole = real(fft(wholedata,axis=0) * conj(fft(wholedata,axis=0)))
        datashape = self.powerspectra_whole.shape[0]
        self.plotfig=figure(fignum)
        if clear: self.plotfig.clear()
        if logy:
            if logx:
                plotcmd = loglog
            else:
                plotcmd = semilogy
        else:
            plotcmd = plot
        plotcmd(fftfreq(datashape,d=self.sample_interval)[0:datashape/2],
                self.powerspectra_whole[0:datashape/2,bolonum],
                linewidth=0.5,color=color)
        xlabel("Frequency (Hz)")

    def broken_powerfit(self, bolonum=0, plbreak=2, doplot=True, logx=True,
            replotspec=True, defaultplot=False, p0in=None, p1in=None, p2in=None, **kwargs):
        if replotspec or (self.powerspectra_whole is None):
            self.powerspec_whole(bolonum=bolonum,logx=True,**kwargs)
        datashape = self.powerspectra_whole.shape[0]
        xfreq = fftfreq(datashape,d=self.sample_interval)[1:datashape/2]
        powerspectra_half = self.powerspectra_whole[1:datashape/2,bolonum]
        p0 = polyfit(log10(xfreq[(xfreq<0.02)]),log10(powerspectra_half[(xfreq<0.02)]),1)
        p1 = polyfit(log10(xfreq[(xfreq<plbreak)*(xfreq>0.02)]),log10(powerspectra_half[(xfreq<plbreak)*(xfreq>0.02)]),1)
        p2 = polyfit(log10(xfreq[xfreq>=plbreak]),log10(powerspectra_half[xfreq>=plbreak]),1)
        # renormalize so that the high-frequency matches the low
        p2nought = p2[1]
        p2[1] = log10(10**p1[1]*(plbreak**(p1[0]-p2[0])))
        def f(x,p):
            return 10**(p[1])*x**(p[0])
        if None not in (p1in,p2in,p0in):
            plot(xfreq[xfreq<0.02],f(xfreq[xfreq<0.02],p0in),color='r')
            plot(xfreq[xfreq<plbreak],f(xfreq[xfreq<plbreak],p1in),color='r')
            plot(xfreq[xfreq>=plbreak],f(xfreq[xfreq>=plbreak],p2in),color='r')
        if doplot: 
            P = plot(xfreq[xfreq<0.02],f(xfreq[xfreq<0.02],p0))[0]
            plot(xfreq[xfreq<plbreak],f(xfreq[xfreq<plbreak],p1),color=P.get_color())
            plot(xfreq[xfreq>=plbreak],f(xfreq[xfreq>=plbreak],p2),color=P.get_color())
        print "Best powerlaw fit: P = 10^%0.3f freq^%0.3f   { freq < %0.2f" % (p1[1],p1[0],plbreak)
        print "                       10^%0.3f freq^%0.3f   { freq >= %0.2f" % (p2[1],p2[0],plbreak)
        print "                       10^%0.3f freq^%0.3f   { freq < %0.2f" % (p0[1],p0[0],0.02)
        print "Reminder: the high-frequency end is forced to meet the low-freqency.  Scale was originally 10^%0.3f" % (p2nought)
        noise_scale = (self.powerspectra_whole[1:datashape/2,0]/f(xfreq,p2)*(xfreq>=2.5))[(9>xfreq)*(xfreq>=2.5)].std()
        print "The Gaussian stddev should be %f and the mean should be 1" % (noise_scale)
        return p0,p1,p2,noise_scale

    def broken_expfit(self,bolonum=0,plbreak=2.5,doplot=True,logx=False,replotspec=True,defaultplot=False,**kwargs):
        """
        Fit two exponentials (one most likely flat) to the power spectrum
        Ignore frequencies < 0.02 Hz, as these are filtered out by the AC sampler
        """
        if replotspec or (self.powerspectra_whole is None):
            self.powerspec_whole(bolonum=bolonum,logx=logx,**kwargs)
        datashape = self.powerspectra_whole.shape[0]
        xfreq = fftfreq(datashape,d=self.sample_interval)[1:datashape/2]
        powerspectra_half = self.powerspectra_whole[1:datashape/2,bolonum]
        p1 = polyfit((xfreq[(xfreq<plbreak)*(xfreq>0.02)]),log10(powerspectra_half[(xfreq<plbreak)*(xfreq>0.02)]),1)
        p2 = polyfit((xfreq[xfreq>=plbreak]),log10(powerspectra_half[xfreq>=plbreak]),1)
        # renormalize so that the high-frequency matches the low
        #p2nought = p2[1]
        #p2[1] = log10(10**p1[1]*(10**(p1[0]*plbreak-p2[0])))
        def f(x,p):
            return 10**(p[0]*x + p[1])
        if defaultplot:
            plot(xfreq[xfreq<plbreak],10**4*(xfreq[xfreq<plbreak])**(-2.0),color='r')
            plot(xfreq[xfreq>=plbreak],2.5e3*(xfreq[xfreq>=plbreak])**(0.0),color='r')
        if doplot: 
            P = plot(xfreq[xfreq<plbreak],f(xfreq[xfreq<plbreak],p1))[0]
            plot(xfreq[xfreq>=plbreak],f(xfreq[xfreq>=plbreak],p2),color=P.get_color())
        print "Best powerlaw fit: P = 10^(%0.3f nu +%0.3f)   { freq < %0.2f" % (p1[0],p1[1],plbreak)
        print "                       10^(%0.3f nu +%0.3f)   { freq >= %0.2f" % (p2[0],p2[1],plbreak)
        noise_scale = (self.powerspectra_whole[1:datashape/2,0]/f(xfreq,p2)*(xfreq>=2.5))[(9>xfreq)*(xfreq>=2.5)].std()
        print "Therefore the Gaussian stddev should be %f and the mean should be 1" % (noise_scale)
        return p1,p2,noise_scale

    def timestream_whole(self,bolonum=0,clear=True,timestream='data',fignum=4, **kwargs):
        if timestream == 'data':
            wholedata = reshape(self.data,[self.data.shape[0]*self.data.shape[1],self.data.shape[2]])
        elif self.tsplot_dict.has_key(timestream):
            data = self.lookup(timestream) #tsplot_dict[timestream]()
            wholedata = reshape(data,[data.shape[0]*data.shape[1],data.shape[2]])
        else:
            raise KeyError("Timestream %s is not valid." % timestream)
        self.plotfig=figure(fignum)
        if clear: self.plotfig.clear()
        if kwargs.has_key('label'):
            label = kwargs.pop('label')
        else:
            label = str(bolonum)
        plot(wholedata[:,bolonum],linewidth=0.5, label=label, **kwargs)


    def ordered_timestreams(self,  ordertype='astro', clear=True,
            colors=['black','purple','blue','cyan','green','orange','red','magenta',],
            astro_order=['ac_bolos','atmo_one','atmo_one_itermedian','atmos_remainder','PCA_astro','astrosignal'],
            atmo_order=['ac_bolos','atmo_one','expmodel','first_sky','PCA_atmo','PCA_astro','noise'],
            dolegend=True, dosubplots=True, fignum=4,
            **kwargs):
        """
        Plot the timestreams in the order they're produced during data reduction

        *ordertype* - 'astro' or 'atmo'
        """
        if ordertype=='astro':
            order = astro_order
        elif ordertype=='atmo':
            order = atmo_order

        figure(fignum)
        clf()

        for ii,tsname in enumerate(order):
            if dosubplots:
                subplot(len(order),1,ii+1)
                title(tsname)
            self.timestream_whole(clear=False, timestream=tsname, label=tsname, color=colors[ii], fignum=fignum, **kwargs)

        if dolegend and not dosubplots:
            L=legend(loc='best')
            L.draggable()
        
    
    def write_ncdf(self):
        if not self.ncfile:
            print "Not writing NCDF file"
            return

#      flags = asarray(ncfile.variables['flags'])
#      bolo_params = asarray(ncfile.variables['bolo_params'])
#      scans_info = asarray(ncfile.variables['scans_info'])
        flags = copy.copy(asarray(self.ncflags))
        scans_info = asarray(self.ncscans)
        bolo_params = asarray(self.ncbolo_params)
        nbolos = self.nbolos
        scanlen = self.scanlen
        nscans = self.nscans
#      self.ngoodbolos = bolo_params[:,0].sum()
        bolo_indices = (self.bolo_indices[newaxis,:] 
                + zeros([self.whscan.shape[0],1]) ).astype('int')
        whscan = (self.whscan[:,newaxis] 
                + zeros([1,self.ngoodbolos])).astype('int')
#      fs= reshape(self.flags,[nscans*scanlen,ngoodbolos])
#      fs2 = zeros([nscans*scanlen,nbolos])
#      fs2[:,self.bolo_indices] = fs
#      flags[self.whscan,:] = fs2
        flags[whscan,bolo_indices] = reshape(self.flags,
                [nscans*scanlen,self.ngoodbolos])
        if flags.min() < 0:
             flags[flags<0] = 0
        #self.ncfile.close()
        #ncfile = netcdf.netcdf_file(self.ncfilename.replace("_ds5","_ds5_flagged"),'w') # NetCDF.NetCDFFile(self.ncfilename,'a')
        #ncfile.variables = self.ncfile.variables
        #ncfile.dimensions = self.ncfile.dimensions
        ncfile = copy.copy(self.ncfile)
        ncfile.filename = ncfile.filename.replace("_ds5","_ds5_flagged")
        ncfile.fp = open(ncfile.filename,'w')
        ncfile.mode = 'w'
        ncfile.variables['flags'].data = flags
        ncfile.createDimension('one',1)
        for key,var in ncfile.variables.items(): 
            if var.shape is ():
                var.dimensions = ('one',)
            if var.__dict__.has_key('file'):
                var.file = (var.file,)
            if var.__dict__.has_key('units'):
                var.units = (var.units,)
            #if type(var.data) in (type(1),type(1.0)):
            #    var.data = [var.data]
            for k,v in var.__dict__.items():
                if k not in ('_shape','_attributes','add_offset','scale_factor'):
                    try:
                        TEMP = v[0]
                    except:
                        print "Failed at %s with value %s and type %s" % (k,v,type(v))
        ncfile.history += "\n Flagged on "+time.ctime()
        ncfile._write()
#      print ncfile.variables['flags'].max()
#      import pdb; pdb.set_trace()
        ncfile.close()

    def redraw(self):
        self.plotscan(self.scannum)

    def _refresh(self):
        if self.flagfig is not None:
            self.flagfig.canvas.draw()
        if self.datafig is not None:
            self.datafig.canvas.draw()
        if self.mapfig is not None:
            self.mapfig.canvas.draw()
        if self.plotfig is not None:
            self.plotfig.canvas.draw()
        self.PCAflag = False
        self.powerspec_plotted = False

    def histograms(self, fignum=8, clear=True, hrange=[0,2], nbins=21, dolegend=True,
            loc='best', ignore_zeros=True, **kwargs):

        self.histfig = figure(fignum)
        if clear: self.histfig.clear()

        if ignore_zeros:
            OKscale = self.scale_coeffs != 0
            OKweight = self.weight_by_bolo != 0
        else:
            OKscale = OKweight = self.scale_coeffs*0+True

        hist(self.scale_coeffs[OKscale],histtype='step',bins=linspace(hrange[0],hrange[1],nbins),label='Scale Coeffs',color='b',linewidth=2)
        hist(self.weight_by_bolo[OKweight],histtype='step',bins=linspace(hrange[0],hrange[1],nbins),label='Weights',color='r',linewidth=2)
        
        if dolegend: legend(loc=loc)

    def unmask_timestream(self, timestream='data'):
        """
        Remove masks for timestream
        """
        if timestream == 'data':
            if hasattr(self.data,'mask'):
                self.data.mask[:] = False
        elif self.tsplot_dict.has_key(timestream):
            data = self.lookup(timestream)
            if hasattr(data,'mask'):
                data.mask[:] = False

    def unmask_all(self):
        for key in self.tscache:
            if hasattr(self.tscache[key],'mask'):
                self.tscache[key].mask[:] = False

    def clearflags(self):
        self.flags[:] = False

        datums=['astrosignal','atmosphere','ac_bolos','atmo_one','noise','scalearr','weight','mapped_astrosignal']
        for d in datums:
            if self.__dict__.has_key(d):
                del self.__dict__[d]
            setattr(self.__class__, d, lazydata(d,flag=False))

        self.unmask_all()

    def doPCA(self,clear=True,timestream='data', fignum=7, plotitem='evects', **kwargs):

        if timestream == 'data':
            self.efuncarr,self.covmat,self.evals,self.evects = efuncs(reshape(self.data,[self.data.shape[0]*self.data.shape[1],self.data.shape[2]]),return_others=True)
        elif self.tsplot_dict.has_key(timestream):
            data = self.lookup(timestream) #tsplot_dict[timestream]()
            self.efuncarr,self.covmat,self.evals,self.evects = efuncs(reshape(data,[data.shape[0]*data.shape[1],data.shape[2]]),return_others=True)
        else:
            raise KeyError("Timestream %s is not valid." % timestream)
        self.PCAfig=figure(fignum)
        if clear: self.PCAfig.clear()
      
        plotitem_dict = {'evects':self.evects,'efuncarr':self.efuncarr,'covmat':self.covmat,'evals':self.evals}
        plottype_dict = {'evects':imshow,'efuncarr':imshow,'covmat':imshow,'evals':plot}
        plotkwargs_dict = {'evects':{},'efuncarr':{'aspect':float(self.efuncarr.shape[1])/self.efuncarr.shape[0]},'covmat':imshow,'evals':plot}
        xlabel_dict = {'evects':'Eigenfunction','efuncarr':'Eigenfunction','covmat':'Bolometer','evals':'Eigenfunction'}
        ylabel_dict = {'evects':'Bolometer','efuncarr':'Time','covmat':'Bolometer','evals':'Correlation Coefficient'}

        plottype_dict[plotitem](plotitem_dict[plotitem], label=plotitem, **dict(plotkwargs_dict[plotitem].items(), kwargs.items()))
        xlabel(xlabel_dict[plotitem])
        ylabel(ylabel_dict[plotitem])
        colorbar()

    def compute_map(self,ts=None,tsname=None,weights=None,showmap=True,**kwargs):
        """
        Create a map from the data and potentially show it
        """
        t0 = time.time()
        if ts is None: 
            if tsname in self.tsplot_dict.keys():
                ts = self.lookup(tsname) #tsplot_dict[tsname]()
            else:
                ts = self.mapped_timestream
        if weights is None: weights = self.weight
        elif not isinstance(weights,numpy.ndarray): weights = numpy.ones(ts.shape)*weights

        self.map = drizzle(self.tstomap,ts,self.map.shape,weights*(True-self.flags))
        print "Computing map took %f seconds" % (time.time() - t0)
        if showmap: self.showmap(**kwargs)

def nantomask(arr):
      mask = (arr != arr)
      return numpy.ma.masked_where(mask,arr)

def masktozero(arr):
      try:
          arr[arr.mask] = 0
      except AttributeError:
          arr[arr!=arr] = 0
      return numpy.array(arr)

def getmask(arr):
      if hasattr(arr,'mask'):
          return arr.mask
      else:
          return numpy.isnan(arr)+(arr==0)

def downsample(myarr,factor):
      factor = int(factor)
      xs,ys = myarr.shape
      crarr = myarr[:xs-(xs % int(factor)),:ys-(ys % int(factor))]
      dsarr = numpy.concatenate([[crarr[i::factor,j::factor]
          for i in xrange(factor)]
          for j in xrange(factor)]).mean(axis=0)
      return dsarr 

def gridmap(x,y,v,downsample_factor=2,smoothpix=3.0,dosmooth=True,xsize=None,ysize=None):
      nx = xrange = numpy.ceil(numpy.max(x)-numpy.min(x))+3
      ny = yrange = numpy.ceil(numpy.max(y)-numpy.min(y))+3
      xax = x-min(x)
      yax = y-min(y)
      if xsize and ysize:
          map = zeros([xsize,ysize])
      else:
          map = zeros([yrange,xrange])
      map[numpy.round(yax),numpy.round(xax)] += v
      map[map!=map] = 0

      if dosmooth:
          dm = smooth(map,kernelwidth=smoothpix)
          #xax,yax = numpy.indices(map.shape)
          #kernel = gaussfitter.twodgaussian([1,nx/2,ny/2,smoothpix],circle=1,rotate=0,vheight=0)(xax,yax)
          #kernelfft = numpy.fft.fft2(kernel)
          #imfft = numpy.fft.fft2(map)
          #dm = numpy.fft.fftshift(numpy.fft.ifft2(kernelfft*imfft).real)
      else:
          dm = map

      return downsample(dm,downsample_factor)

def _hdr_string_to_card(str):
      name = str[:7].strip()
      val  = str[9:31].strip()
      try:
          val = float(val)
      except:
          pass
      comment = str[31:].strip()
      if name == 'END':
          return
      else:
          return pyfits.Card(name,val,comment)

def _hdr_string_list_to_cardlist(strlist):
      cardlist = [_hdr_string_to_card(s) for s in strlist]
      while None in cardlist: cardlist.remove(None)
      return pyfits.CardList(cardlist)

def itermedian(ts,scale=0.8,niter=5,axis=2):
      mts = ts.copy()
      for i in xrange(niter):
          mts = mts - numpy.median(mts,axis=2)[:,:,newaxis]*scale
      return mts

def boloexp(xax,height,amp,fjac=None): 
    # set time constant of exponential:
    tc=10.4 #in seconds
    expterm = exp(-xax/(tc*50.))	#data taken at 50 Hz
    model = amp * expterm + height
    return model

def boloexp_resids(xax,data):
    def f(p,fjac=None):
        model = boloexp(xax,p[0],p[1])
        return [0,(data-model)]
    return f

def expsub_line(ts, sampleinterval=0.04, quiet=True, **kwargs):
    xax = arange(ts.shape[0]) * sampleinterval/0.02
    mp = mpfit.mpfit(boloexp_resids(xax,ts),[0,median(ts)],quiet=quiet,**kwargs)
    model = boloexp(xax,*mp.params)
    return ts-model

def exponent_sub(arr, **kwargs):
    nscans,ntime,nbolos = arr.shape
    subbedarr = np.zeros(arr.shape)
    for ii in xrange(nbolos):
        for jj in xrange(nscans):
            if hasattr(arr,'mask'):
                # make sure there's enough data to fit
                if arr.mask[jj,:,ii].sum() < ntime - 2:
                    subbedarr[jj,:,ii] = expsub_line(arr[jj,:,ii], **kwargs)
            else:
                subbedarr[jj,:,ii] = expsub_line(arr[jj,:,ii], **kwargs)

    return subbedarr


if __name__ == "__main__":

      #from pylab import *
      #from agpy import pyflagger
      import sys
      import optparse
      import pdb

      parser=optparse.OptionParser()
      parser.add_option("--timestreams","-t",help="Plot all timestreams scaled and unscaled as a diagnostic tool.  Default False",action="store_true",dest="timestreams",default=False)
      parser.add_option("--plotscan","-p",help="Scan number to plot at the start.  Default to 0.  Set to -1 to turn off.",default=0)
      parser.add_option("--interactive","-i",help="Start in ipython mode?  Default True",action='store_false',default=True)
      #parser.add_option("--debug","-d",help="Debug in ipython mode?",action='store_true',default=False)
      parser.add_option("--no_interactive","-n",help="Turn off ipython (interactive) mode",action='store_false',dest="interactive")
      parser.add_option("--close","-c",help="Close after plotting?  Useful if interactive=False",action='store_true',default=False)
      parser.add_option("--noflag",help="Disable flagging?",action='store_true',default=False)
      parser.add_option("--compute_expfit",help="Do exponential fits to the powerspectra for all bolometers?",action='store_true',default=False)
      parser.add_option("--compute_powerfit",help="Do power-law fits to the powerspectra for all bolometers?",action='store_true',default=False)
      parser.set_usage("%prog savename.sav [options]")
      parser.set_description(
      """
      PYFLAGGER 
      """)

      options,args = parser.parse_args()

      flag = not options.noflag
          
      plotnum = options.plotscan

      flist = []
      for ii,filename in enumerate(args):
          prefix = filename.replace("_preiter.sav","").replace("_postiter.sav","").replace("_save_iter0.sav","")

          try:
              f=Flagger(filename, flag=flag)
          except:
              pdb.post_mortem()
          if options.interactive: exec("f%i = f" % (ii) )

          if plotnum >= 0: f.plotscan(plotnum)

          if options.timestreams:
              figure(6)
              clf()
              title("ac_bolos scaled")
              for ii in range(f.nbolos):
                  f.timestream_whole(ii,timestream='acbolos',clear=False,fignum=6)
              savefig('%s_acbolos_scaled.png' % prefix)

              figure(5)
              clf()
              title("ac_bolos unscaled")
              for ii in range(f.nbolos):
                  f.timestream_whole(ii,timestream='acbolos_noscale',clear=False,fignum=5)
              savefig('%s_acbolos_noscale.png' % prefix)

          if options.compute_powerfit:
              ppars0 = zeros([2,f.nbolos])
              ppars1 = zeros([2,f.nbolos])
              ppars2 = zeros([2,f.nbolos])
              gstd = zeros(f.nbolos)
              for ii in range(f.nbolos):
                  try:
                      ppars0[:,ii],ppars1[:,ii], ppars2[:,ii], gstd[ii] = f.broken_powerfit(ii,timestream='ac_bolos')
                  except ValueError:
                      print "Bolo %i is probably flagged out." % ii
              print "Median(gstd) = %f" % median(gstd)
              print "Median(ppars0[0,:]) = %f" % median(ppars0[0,:])
              print "Median(ppars0[1,:]) = %f" % median(ppars0[1,:])
              print "Median(ppars1[0,:]) = %f" % median(ppars1[0,:])
              print "Median(ppars1[1,:]) = %f" % median(ppars1[1,:])
              print "Median(ppars2[0,:]) = %f" % median(ppars2[0,:])
              print "Median(ppars2[1,:]) = %f" % median(ppars2[1,:])
              print "WARNING! This likely only applies to astrosignal, NOT ac_bolos!"

          if options.compute_expfit:
              ppars1 = zeros([2,f.nbolos])
              ppars2 = zeros([2,f.nbolos])
              gstd = zeros(f.nbolos)
              for ii in range(f.nbolos):
                  try:
                      ppars1[:,ii], ppars2[:,ii], gstd[ii] = f.broken_expfit(ii,timestream='ac_bolos')
                  except ValueError:
                      print "Bolo %i is probably flagged out." % ii
              print "Median(gstd) = %f" % median(gstd)
              print "Median(ppars1[0,:]) = %f" % median(ppars1[0,:])
              print "Median(ppars1[1,:]) = %f" % median(ppars1[1,:])
              print "Median(ppars2[0,:]) = %f" % median(ppars2[0,:])
              print "Median(ppars2[1,:]) = %f" % median(ppars2[1,:])
              print "WARNING! This likely only applies to astrosignal, NOT ac_bolos!"

          if options.close or not options.interactive:
              f.close()
          else:
              flist.append(f)


      if options.interactive:
          try:
              from IPython.Shell import IPShellEmbed
              ipshell = IPShellEmbed([])
              ipshell() # this call anywhere in your program will start IPython
          except ImportError:
              from IPython import embed
              embed()

