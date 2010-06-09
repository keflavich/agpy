
import math

import pylab
from pylab import *
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
import copy
import idlsave
import gaussfitter
from PCA_tools import *

matplotlib.rcParams['image.origin']='lower'
matplotlib.rcParams['image.interpolation']='nearest'
matplotlib.rcParams['image.aspect']=1
matplotlib.defaultParams['image.origin']='lower'
matplotlib.defaultParams['image.interpolation']='nearest'
matplotlib.defaultParams['image.aspect']=1
if matplotlib.rcParams['text.usetex']:
  texOn = True
# matplotlib.rcParams['text.usetex']=False
# matplotlib.defaultParams['text.usetex']=False

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
 
  Map Key Commands:
    c - toggle current scan
    . - show point in timestream
    click - show point in timestream
    r - redraw

  """

  def __init__(self, filename, debug=False, **kwargs):
      self.filename = filename
      self.debug = debug
      if filename[-4:] == 'fits':
          self._loadfits(filename,**kwargs)
      elif filename[-3:] == 'sav':
          self._loadsav(filename,**kwargs)

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
      self.PCA = False

      self.showmap(vmax=vmax)
      self.dcon()
  
  def _loadsav(self,savfile,**kwargs):
      print "Beginning IDLsave file read."
      t0 = time.time()
      self.sav = idlsave.read(savfile)
      print "Completed IDLsave file read in %f seconds." % (time.time() - t0)

      self.ncfilename = savfile
      self.tsfile = None

      self.ncscans = self.sav.variables['bgps']['scans_info'][0]
      if len(self.ncscans.shape) == 1: self.ncscans.shape = [1,2]
      self.scanlengths = self.ncscans[:,1]+1-self.ncscans[:,0]
      self.scanlen = numpy.max(self.scanlengths)
      self.ncflags = self.sav.variables['bgps']['flags'][0] 
      self.timelen = self.ncflags.shape[0]
      self.nbolos = self.ncflags.shape[1]
      self.nscans = self.ncscans.shape[0]
      self.ncbolo_params = self.sav.variables['bgps']['bolo_params'][0]
      self.ncbolo_indices = self.sav.variables['bgps']['bolo_indices'][0]
      #self.bolo_indices = asarray(nonzero(self.ncbolo_params[:,0].ravel())).ravel()
      self.bolo_indices = self.ncbolo_indices
      self.ngoodbolos = self.bolo_indices.shape[0]
      self.whscan = asarray([arange(self.scanlen)+i for i,j in self.ncscans]).ravel()
      self.scanstarts = arange(self.nscans)*self.scanlen
      self.whempty = concatenate([arange(i+j,i+self.scanlen) for i,j in zip(self.scanstarts,self.scanlengths) ]).ravel()
      self.whscan[self.whempty] = 0

      self.datashape = [self.nscans,self.scanlen,self.ngoodbolos]

      self.astrosignal = self.sav.variables['bgps']['astrosignal'][0][self.whscan,:].astype('float')
      self.atmosphere  = self.sav.variables['bgps']['atmosphere'][0][self.whscan,:].astype('float')
      self.raw         = self.sav.variables['bgps']['raw'][0][self.whscan,:].astype('float')
      self.ac_bolos    = self.sav.variables['bgps']['ac_bolos'][0][self.whscan,:].astype('float')
      self.atmo_one    = self.sav.variables['bgps']['atmo_one'][0][self.whscan,:].astype('float')
      self.dc_bolos    = self.sav.variables['bgps']['dc_bolos'][0][self.whscan,:].astype('float')
      self.noise       = self.sav.variables['bgps']['noise'][0][self.whscan,:].astype('float')
      self.scalearr    = self.sav.variables['bgps']['scalearr'][0][self.whscan,:].astype('float')
      self.weight      = self.sav.variables['bgps']['weight'][0][self.whscan,:].astype('float')
      self.flags       = self.sav.variables['bgps']['flags'][0][self.whscan,:]
      self.flags.shape = self.datashape

      try:
          self.mapped_astrosignal = self.sav.variables['bgps']['mapped_astrosignal'][0][self.whscan,:].astype('float')
      except:
          self.mapped_astrosignal = copy.copy(self.astrosignal)

      datums=['astrosignal','atmosphere','raw','ac_bolos','atmo_one','dc_bolos','noise','scalearr','weight','mapped_astrosignal']
      for d in datums:
          self.__dict__[d][self.whempty,:] = NaN
          self.__dict__[d].shape = self.datashape
          self.__dict__[d] = nantomask(self.__dict__[d])
          try:
              self.__dict__[d].mask[self.flags > 0] = True
          except TypeError:
              self.__dict__[d].mask = (self.flags > 0)
      #import pdb; pdb.set_trace()

      if list(self.ac_bolos.shape) != self.datashape:
          import pdb; pdb.set_trace()
      self.data = self.ac_bolos

      self.ncfile = None
      self.flagfn = savfile.replace("sav","_flags.fits")

      self.map      = nantomask( self.sav.variables['mapstr']['astromap'][0] )
      self.tstomap  = reshape( self.sav.variables['mapstr']['ts'][0][self.whscan,:] , self.datashape )

      self.mapped_timestream = self.atmo_one - self.atmosphere + self.astrosignal

      if self.map.sum() == 0:
          self.map  = nantomask( self.sav.variables['mapstr']['rawmap'][0] )

      self.header = pyfits.Header(_hdr_string_list_to_cardlist( self.sav.variables['mapstr']['hdr'][0] ))

      self._initialize_vars(**kwargs)

      self.tsplot = 'default'
      self.set_tsplot(**kwargs)

  
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
      """
      if tsplot is not None:
          self.tsplot=tsplot
      if self.tsplot == 'astrosignal' and self.astrosignal.sum() != 0:
          self.data = self.astrosignal
      elif self.astrosignal.sum() == 0:
          print "Astrosignal is zero."
      elif self.tsplot == 'dcbolos':
          self.data = self.dc_bolos*self.scalearr
      elif self.tsplot == 'acbolos_noscale' or self.tsplot == 'ac_bolos_noscale':
          self.data = self.ac_bolos
      elif self.tsplot == 'atmo_one':
          self.data = self.atmo_one
      elif self.tsplot == 'acbolos' or self.tsplot=='ac_bolos':
          self.data = self.ac_bolos*self.scalearr
      elif self.tsplot == 'atmosphere':
          self.data = self.atmosphere
      elif self.tsplot=='skysub_noscale':
          self.data = self.atmo_one - self.atmosphere
      elif self.tsplot=='residual':
          self.data = self.atmo_one - self.atmosphere - self.noise
      elif self.tsplot=='skysub' or self.tsplot=='default':
          self.data = self.atmo_one - self.atmosphere + self.astrosignal
      elif self.tsplot=='last_astrosignal':
          self.data = self.atmo_one - self.atmosphere - self.noise + self.astrosignal
      elif self.tsplot=='default_noscale':
          self.data = self.ac_bolos - self.atmo_one - self.atmosphere
      elif self.tsplot=='scale':
          self.data = self.scalearr
      elif self.tsplot=='weight':
          self.data = self.weight
      elif self.tsplot=='raw':
          self.data = self.raw
      elif self.tsplot=='rawscaled':
          self.data = self.raw * self.scalearr
      elif self.tsplot=='noise':
          self.data = self.noise
      elif self.tsplot=='mapped_astrosignal':
          self.data = self.mapped_astrosignal
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


  def showmap(self,colormap=cm.spectral,vmin=None,vmax=None):
    self.mapfig=figure(0); clf(); 
    if vmax is None:
        vmax = self.map.mean()+7*self.map.std()
    elif vmax=='max':
        vmax = self.map.max()
    if vmin is None:
        vmin = self.map.mean()-2*self.map.std()
    elif vmin=='min':
        vmin = self.map.min()
    self.mapim = pylab.imshow(self.map,
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
      try:
          disconnect(self.gfit_connection)
      except:
          pass
      self.gfit_connection = connect('button_press_event',self.gfit_bolomap)
      pylab.colorbar()
      return self.bolommap

  def gfit_map(self,event,map,ax=None):
      if event.xdata and event.ydata:
          gf = gaussfitter.gaussfit(masktozero(map),
                  params=   [0,0,event.xdata,event.ydata,0,0,0],
                  usemoment=[1,1,          0,          0,1,1,1],
                  rotate=1)
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

  def footmovie(self,y1,y2,movie=False,moviedir='scanmovie/',logscale=False,smooth=True):
      y1 = numpy.round(y1)
      y2 = numpy.round(y2)
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
      mapcube = array([ gridmap(x[ii,:],y[ii,:],flagvals[ii,:],downsample_factor=2,xsize=72,ysize=72,smooth=smooth)
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
        xlabel('Bolometer number'); ylabel('Time (0.1s)')
        self.flagim = pylab.imshow(self.flags[scannum,:,:],interpolation='nearest',
                origin='lower',aspect=self.aspect)
        self.flagcb = pylab.colorbar()
        self.datafig = figure(fignum,figsize=[16,12]);clf();
        if texOn:
          self.datatitle = pylab.title("Scan "+str(self.scannum)+" in "+self.ncfilename.replace("_","\\_"));
        else:
          self.datatitle = pylab.title("Scan "+str(self.scannum)+" in "+self.ncfilename);
        xlabel('Bolometer number'); ylabel('Time (0.1s)')
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
      if button=='b' and self.PCA:
          x=round(x)
          PCAsub = self.data[self.scannum,:,:] - outer(self.plane[:,x],ones(self.plane.shape[1]))
          self.PCA = False
          self.plotscan(self.scannum,data=PCAsub)
      elif button=='b' and not self.PCA:
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

  def plot_column(self,tsx):
      self.bolofig=figure(4)
      self.bolofig.clf()
      pylab.plot(self.plane[:,numpy.round(tsx)])

  def plot_line(self,tsy):
      self.bolofig=figure(4)
      self.bolofig.clf()
      if self.PCA:
          pylab.plot(self.efuncarr[numpy.round(tsy),:])
      else:
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
      tempdata = self.data
      if self.tsfile:
          self.tsfile[0].data = asarray(self.flags,dtype='int')
          self.tsfile.writeto(self.flagfn,clobber=True)
          self.tsfile[0].data = tempdata
      elif self.header:
          self.flagfits = pyfits.PrimaryHDU(asarray(self.flags,dtype='int'),self.header)
          self.flagfits.writeto(self.flagfn,clobber=True)
          print "flags_to_ncdf,'%s','%s'" % (self.flagfn,self.filename)

  def mapclick(self,event):
      if event.xdata == None:
          return
      clickX = round(event.xdata)
      clickY = round(event.ydata)
      if event.button == 1:
          self.tsarrow(clickX,clickY)
      elif event.button == 2:
          self.find_all_points(clickX,clickY)

  def mapkeypress(self,event):
      if event.inaxes is None: return
      elif event.key == 'c':
          self.toggle_currentscan()
      elif event.key == 'G':
          self.gfit_map(event,self.map,self.mapim.axes)
          self.mapfig.canvas.draw()
      elif event.key == '.':
          if event.xdata == None:
              return
          clickX = round(event.xdata)
          clickY = round(event.ydata)
          self.tsarrow(clickX,clickY)
      elif event.key == "r":
          self.showmap()

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
      elif event.key == 'P': # PCA
          self.plotscan(self.scannum,data=efuncs(self.plane),flag=False,logscale=True)
          self.PCA = True
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
      elif event.key == 'S':
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
          print "Value at %i,%i: %f  Flagged=%i" % (x,y,vpt,fpt)
      if set_lastkey: 
          self._lastkey = event.key

  def find_all_points(self,x,y):
      mappoint = y * self.map.shape[1] + x
      self.timepoints =  nonzero(self.tstomap == mappoint)
      wtavg = (self.mapped_timestream[self.timepoints]*self.weight[self.timepoints]).sum() / self.weight[self.timepoints].sum()
      uwtavg = self.mapped_timestream[self.timepoints].mean()
      print "Map value: %f   Weighted average: %f   Unweighted Average: %f" % (self.map[y,x],wtavg,uwtavg)
      print "scan,bolo,time: %12s%12s%12s%12s%12s" % ('mapped','astro','flags','weight','scale')
      for ii,jj,kk in transpose(self.timepoints):
          print "%4i,%4i,%4i: %12f%12f%12f%12f%12f" % (ii,kk,jj,
                  self.mapped_timestream[ii,jj,kk],
                  #self.mapped_astrosignal[ii,jj,kk],
                  self.astrosignal[ii,jj,kk],
                  self.flags[ii,jj,kk],
                  self.weight[ii,jj,kk],
                  self.scalearr[ii,jj,kk])


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
    if tb.mode=='' and not self.PCA:
        self.flag_box(self.x1,self.y1,self.x2,self.y2,event.button)
#        if abs(self.x2-self.x1) > 1 or abs(self.y2-self.y1) > 1:
#        else:
#            self.flagpoint(self.x1,self.y1,event.button)

  def mouse_down_event(self, event):
    if event.inaxes is None: return
    self.mouse_up = False
    self.x1 = event.xdata
    self.y1 = event.ydata

  
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

  def _refresh(self):
      if self.flagfig is not None:
          self.flagfig.canvas.draw()
      if self.datafig is not None:
          self.datafig.canvas.draw()
      if self.mapfig is not None:
          self.mapfig.canvas.draw()
      if self.plotfig is not None:
          self.plotfig.canvas.draw()
      self.PCA = False

def nantomask(arr):
    mask = (arr != arr)
    return numpy.ma.masked_where(mask,arr)

def masktozero(arr):
    try:
        arr[arr.mask] = 0
    except AttributeError:
        arr[arr!=arr] = 0
    return numpy.array(arr)

def downsample(myarr,factor):
    factor = int(factor)
    xs,ys = myarr.shape
    crarr = myarr[:xs-(xs % int(factor)),:ys-(ys % int(factor))]
    dsarr = numpy.concatenate([[crarr[i::factor,j::factor]
        for i in xrange(factor)]
        for j in xrange(factor)]).mean(axis=0)
    return dsarr 

def gridmap(x,y,v,downsample_factor=2,smoothpix=3.0,smooth=True,xsize=None,ysize=None):
    nx = xrange = numpy.ceil(numpy.max(x)-numpy.min(x))+3
    ny = yrange = numpy.ceil(numpy.max(y)-numpy.min(y))+3
    xax = x-min(x)
    yax = y-min(y)
    if xsize and ysize:
        map = zeros([xsize,ysize])
    else:
        map = zeros([yrange,xrange])
    map[numpy.round(yax),numpy.round(xax)] += v

    if smooth:
        xax,yax = numpy.indices(map.shape)
        kernel = gaussfitter.twodgaussian([1,nx/2,ny/2,smoothpix],circle=1,rotate=0,vheight=0)(xax,yax)
        kernelfft = numpy.fft.fft2(kernel)
        imfft = numpy.fft.fft2(map)
        dm = numpy.fft.fftshift(numpy.fft.ifft2(kernelfft*imfft).real)
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
    cardlist.remove(None)
    return pyfits.CardList(cardlist)

