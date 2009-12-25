
import math

import pylab
from pylab import *
import matplotlib
import pyfits
from numpy import isnan
from mad import MAD,nanmedian
#from matplotlib import patches
from matplotlib.patches import Rectangle,FancyArrow
from matplotlib.lines import Line2D
from matplotlib.widgets import Cursor, MultiCursor
import matplotlib.cm as cm
#from Scientific.IO import NetCDF
from scipy.io import netcdf
import time
import re
import os
import copy


matplotlib.defaultParams['image.origin']='lower'
matplotlib.defaultParams['image.interpolation']='nearest'
matplotlib.defaultParams['image.aspect']=.1

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

  def __init__(self, filename, ncfilename='', flagfile='', mapnum='', axis=None, vmax=None):
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
    self.map[isnan(self.map)] = 0
    self.tstomapfn = self.pathprefix+self.fileprefix+"_tstomap.fits"
    self.tstomapfile = pyfits.open(self.tstomapfn)
    self.tstomap = self.tstomapfile[0].data

#    self.outfile = open(self.pathprefix+"log_"+self.fileprefix+"_flags.log",'a')
    self.data = self.tsfile[0].data
    self.aspect = float(self.data.shape[2])/float(self.data.shape[1])
    self.maxscan = self.data.shape[0]
    self.flagfn = self.pathprefix+self.fileprefix+"_flags.fits"
#    if os.path.exists(self.flagfn):
#        self.flagfile = pyfits.open(self.flagfn)
#        self.flags = self.flagfile[0].data
#        if self.flags.shape != self.data.shape:
#            print "Flags / data shape are different.",self.flags.shape,self.data.shape
#    else:
#        self.flagfile = copy.copy(self.tsfile)
#        self.flags = zeros(self.data.shape,dtype='int')
    print "There are %i scans" % (self.data.shape[0])
#    print >>self.outfile,"Started a new session at "+time.asctime()
    self.reset()
    self.counter = 0
    self.mouse_up = False
    self.connected = 0
    #self.renderer = matplotlib.backends.backend_agg.RendererAgg
    self.rectangles=[[] for i in xrange(self.maxscan)]
    self.lines=[[] for i in xrange(self.maxscan)]
    self.arrows=[]
    self.maparrows=[]
    self.md = 0
    self.mu = 0
    self.key = 0
    self.scannum = 0
    self.fignum = 1
    self.open = 1
    self.currentscan = 0

    self.showmap(vmax=vmax)

  def showmap(self,colormap=cm.spectral,vmax=None):
    self.fig0=figure(0); clf(); 
    if not vmax:
        vmax = self.map.mean()+7*self.map.std()
    elif vmax=='max':
        vmax = self.map.max()
    imshow(self.map,vmin=self.map.mean()-2*self.map.std(),
            vmax=vmax,interpolation='nearest',
            cmap=colormap); 
    colorbar()
    try:
        disconnect(self.MtoT)
        disconnect(self.MtoTkey)
    except:
        pass
    self.MtoT = connect('button_press_event',self.mapclick)
    self.MtoTkey = connect('key_press_event',self.mapkeypress)
    self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)
    self.connections=[self.MtoT]

  def __call__(self, event):
    """ this section is never invoked """
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      print event.xdata, event.x
#      tb = get_current_fig_manager().toolbar
#      self.flagpoint(clickY,clickX,event.button)
#      print clickX,clickY
#      print self.axis
#      print event.inaxes
#      if ((self.axis is None) or (self.axis==event.inaxes)) and tb.mode=='':
#        self.flagpoint(clickY,clickX,event.button)


  def footprint(self,tsx,tsy):
    mappoints = asarray(self.tstomap[self.scannum,tsy,:])

#    for a in self.maparrows:
#        a.set_visible(False)
#    for a in self.maparrows:
#        self.maparrows.remove(a)

    x,y = mappoints / self.map.shape[1],mappoints % self.map.shape[1]

    try:
#        self.fp1.set_visible(False)
        self.fp2[0].set_visible(False)
        self.fp3[0].set_visible(False)
#        draw()
    except:
        pass
#        print "Could not set_visible(False)"
#        draw()

    figure(0)
    myaxis = self.fig0.axes[0].axis()
    self.fp3 = plot(y,x,'ro')
#    self.fp1 = plot(y,x,'b+')
    self.fp2 = plot(y,x,'wx')
    self.fig0.axes[0].axis(myaxis)
#    self.fp1 = scatter(y,x,30,'blue','+',edgecolor='blue')
#    self.fp2 = scatter(y,x,30,'red','x',edgecolor='red')
#    self.fig0.axes[0].images[0].cmap =
#    gray()
#    for mp in mappoints:
#        a1 = arrow(y+2,x+2,-4,-4,head_width=0,facecolor='red',
#           edgecolor='red',length_includes_head=True,head_starts_at_zero=False)
#        a2 = arrow(y-2,x+2,4,-4,head_width=0,facecolor='red',
#           edgecolor='red',length_includes_head=True,head_starts_at_zero=False)
#        a1.set_visible(True)
#        a2.set_visible(True)
#        self.maparrows.append(a1)
#        self.maparrows.append(a2)
    draw()

 
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


  def plotscan(self, scannum, fignum=1, button=1):
    if self.connected:
        self.dcon()
    self.connected = 1
    self.plane = self.data[scannum,:,:] * (self.flags[scannum,:,:]==0)
    self.scannum = scannum
    self.fignum = fignum
    self.flagfig = figure(fignum+1)
    clf()
    #subplot(122)
    title("Flags for Scan "+str(self.scannum)+" in "+self.ncfilename);
    xlabel('Bolometer number'); ylabel('Time (.02s)')
    imshow(self.flags[scannum,:,:],interpolation='nearest',
            origin='lower',aspect=self.aspect)
    colorbar()
    self.datafig = figure(fignum);clf(); #subplot(121)
    title("Scan "+str(self.scannum)+" in "+self.ncfilename);
    xlabel('Bolometer number'); ylabel('Time (.02s)')
    imshow(self.plane,interpolation='nearest',
            origin='lower',aspect=self.aspect)
    colorbar()
    self.showrects()
    self.showlines()
    self.cursor = Cursor(gca(),useblit=True,color='black',linewidth=1)
    draw()
    self.md=connect('button_press_event',self.mouse_down_event)
    self.mu=connect('button_release_event',self.mouse_up_event)
    self.key = connect('key_press_event',self.keypress)
    self.connections.append(self.md)
    self.connections.append(self.mu)
    self.connections.append(self.key)

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
          draw()
 
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
          draw()
#          print x,y,w,h,p
#          print >>self.outfile,\
#                  "flag_manual,'%s',bolorange=[%i,%i],timerange=[%i,%i],scanrange=%i" \
#                  % (self.ncfilename,x1,x2,y1,y2,self.scannum)
      elif button==2:
         # this won't work right; I need some way to make it undo-able  
         # <--- I don't know if that's true any more (11/10/08)
          self.flags[scannum,yrange[0]:yrange[1],xrange[0]:xrange[1]] *= (
                  -1*(self.flags[scannum,yrange[0]:yrange[1],xrange[0]:xrange[1]] > 0) )
          p = matplotlib.patches.Rectangle(xy=(x2,y2), width=w, height=h,
                  facecolor='blue',transform=gca().transData)
          gca().add_patch(p)
          p.set_visible(True)
          p.set_alpha(.5)
          self.rectangles[self.scannum].append(p)
          draw()
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
          draw()
      elif button==3:
          self.maparrow(x2i,y2i)
 
  def flag_bolo(self,x,button):
      if button=='b':
#          print x,round(x),button,self.data.shape
          x=round(x)
          h=self.data.shape[1]
          self.flags[self.scannum,0:h,x] += 1
          p = matplotlib.lines.Line2D([x,x],[0,h],\
                  color='black',transform=gca().transData)
          gca().add_line(p)
          p.set_visible(True)
          self.lines[self.scannum].append(p)
          draw()
#          print >>self.outfile,\
#                  "flag_manual,'%s',bolorange=[%i],scanrange=%i" \
#                  % (self.ncfilename,x,self.scannum)

  def flag_time(self,y,button):
      if button=='t':
          y=round(y)
          w=self.data.shape[2]
          self.flags[self.scannum,y,0:w] += 1
          p = matplotlib.lines.Line2D([0,w],[y,y],\
                  color='black',transform=gca().transData)
          gca().add_line(p)
          p.set_visible(True)
          self.lines[self.scannum].append(p)
          draw()
#          print >>self.outfile,\
#                  "flag_manual,'%s',timerange=[%i,%i],scanrange=%i" \
#                  % (self.ncfilename,0,w,self.scannum)

  def unflag_time(self,y,button):
      if button=='T':
          y=round(y)
          w=self.data.shape[2]
          for l in self.lines[self.scannum]:
              if l._y[0] == y and l._y[1] == y:
                  self.flags[self.scannum,y,0:w] -= 1
                  l.set_visible(False)
                  self.lines[self.scannum].remove(l)
                  draw()
          if self.flags[self.scannum,y,0:w].max() > 0:
              arr = self.flags[self.scannum,y,0:w]
              arr[arr>0] -= 1

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
                  draw()
          if self.flags[self.scannum,0:h,x].max() > 0:
              arr = self.flags[self.scannum,0:h,x]
              arr[arr>0] -= 1


  def dcon(self):
      self.connected = 0
      disconnect(self.md)
      disconnect(self.mu)
      disconnect(self.key)
      disconnect(self.MtoT)
      for i in self.connections:
          self.fig0.canvas.mpl_disconnect(i)
          try:
              self.datafig.canvas.mpl_disconnect(i)
              self.flagfig.canvas.mpl_disconnect(i)
          except:
              continue

  def close(self,write=True):
      """ close the ncdf file and the graphics windows
      and flush everything to file"""
      if self.open == 1:
           self.open = 0
           self.dcon()
           if write:
               self.write_ncdf()
#           self.outfile.close()
     #      self.writeflags()
           figure(self.fignum-1); clf()
           figure(self.fignum);   clf()
           figure(self.fignum+1); clf()

  def writeflags(self):
      tempdata = self.data
      self.tsfile[0].data = self.flags
      self.tsfile.writeto(self.flagfn,clobber=True)
      self.tsfile[0].data = tempdata

  def mapclick(self,event):
      if event.xdata == None:
          return
      clickX = round(event.xdata)
      clickY = round(event.ydata)
      self.tsarrow(clickX,clickY)

  def mapkeypress(self,event):
      if event.inaxes is None: return
      elif event.key == 'c':
          self.toggle_currentscan()
      elif event.key == '.':
          if event.xdata == None:
              return
          clickX = round(event.xdata)
          clickY = round(event.ydata)
          self.tsarrow(clickX,clickY)
      elif event.key == "r":
          self.showmap()

  def keypress(self,event):
      if event.inaxes is None: return
      elif event.key == 'n':
          if self.scannum < self.maxscan-1:
              self.plotscan(self.scannum+1)
          else:
              print "At last scan, can't go further"
      elif event.key == 'p':
          if self.scannum > 0:
              self.plotscan(self.scannum-1)
          else:
              print "At first scan, can't go further back"
      elif event.key == 'q':
          self.close()
      elif event.key == 'Q':
          self.close(write=False)
      elif event.key == '.':
          self.maparrow(round(event.xdata),round(event.ydata))
      elif event.key == 'f':
          self.footprint(round(event.xdata),round(event.ydata))
      elif event.key == 'R': # reverse order of boxes
          self.rectangles[self.scannum].reverse()
      elif event.key == 'r': # redraw
        #figure(self.fignum+1); clf()
        #subplot(122)
        #imshow(self.flags[self.scannum,:,:],\
        #interpolation='nearest',origin='lower',aspect=.1)
        #colorbar()
        #figure(self.fignum);clf(); #subplot(121)
        #self.plane = self.data[self.scannum,:,:]\
        #        * (self.flags[self.scannum,:,:]==0)
        #imshow(self.plane,interpolation='nearest',origin='lower',aspect=.1)
        #colorbar()
        self.plotscan(self.scannum)
        #self.showrects()
        #self.showlines()
        #draw()
      elif event.key == 'd':
          self.flag_box(self.x1,self.y1,self.x2,self.y2,'d')
      elif event.key == 't':
          self.flag_time(event.ydata,event.key)
      elif event.key == 's' or event.key == 'w':
          self.flags[self.scannum,:,:] += 1
      elif event.key == 'S':
          self.flags[self.scannum,:,:] -= (self.flags[self.scannum,:,:] > 0)
      elif event.key == 'b':
          self.flag_bolo(event.xdata,event.key)
      elif event.key == 'T':
          self.unflag_time(event.ydata,event.key)
      elif event.key == 'B':
          self.unflag_bolo(event.xdata,event.key)
      elif event.key == 'c':
          self.toggle_currentscan()
      elif event.key == 'v':
          x,y = round(event.xdata),round(event.ydata)
          vpt = self.data[self.scannum,y,x]
          fpt = self.flags[self.scannum,y,x]
          print "Value at %i,%i: %f  Flagged=%i" % (x,y,vpt,fpt)

  def tsarrow(self,x,y):
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
#              print "i shape: ",i.shape, " matchpts ",matchpts
              i = int(i)
              t,b = self.timepoints[1][i],self.timepoints[2][i]
#              print "T,b,i  ",t,b,i
#              print "Does t = []?",t == []
#              print "Is t >= 0?",t >= 0
#              arrow = FancyArrow(t-5,b-5,5,5)
#              self.datafig.axes[0].add_patch(arrow)
              figure(self.fignum)
              self.datafig.sca(self.datafig.axes[0])
              #arrow = self.datafig.axes[0].arrow(t-5,b-5,5,5)
              a1 = arrow(b-3,t-3,6,6,head_width=0,facecolor='black')
              a2 = arrow(b-3,t+3,6,-6,head_width=0,facecolor='black')
              a1.set_visible(True)
              a2.set_visible(True)
#              print a,t,b
              self.arrows.append(a1)
              self.arrows.append(a2)
#          draw()
          self.datafig.canvas.draw()
          self.datafig.canvas.draw()
#          self.datafig

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
    a1 = arrow(y+2,x+2,-4,-4,head_width=0,facecolor='black',
            length_includes_head=True,head_starts_at_zero=False)
    a2 = arrow(y-2,x+2,4,-4,head_width=0,facecolor='black',
            length_includes_head=True,head_starts_at_zero=False)
    a1.set_visible(True)
    a2.set_visible(True)
    self.maparrows.append(a1)
    self.maparrows.append(a2)
    draw()

  def toggle_currentscan(self):
      if self.currentscan == 0:
          xarr = self.tstomap[self.scannum,:,:] / self.map.shape[1]
          yarr = self.tstomap[self.scannum,:,:] % self.map.shape[1]
          x0,x1 = xarr.min(),xarr.max()
          y0,y1 = yarr.min(),yarr.max()
          self.fig0.axes[0].axis([y0,y1,x0,x1])
          self.currentscan = 1
          self.mapcursor=Cursor(gca(),useblit=True,color='black',linewidth=1)
      elif self.currentscan == 1:
          self.fig0.axes[0].axis([0,self.map.shape[1],0,self.map.shape[0]])
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
    if tb.mode=='':
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

#      flags = asarray(ncfile.variables['flags'])
#      bolo_params = asarray(ncfile.variables['bolo_params'])
#      scans_info = asarray(ncfile.variables['scans_info'])
      flags = asarray(self.ncflags)
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
      self.ncfile.close()
      ncfile = netcdf.netcdf_file(self.ncfilename,'w') # NetCDF.NetCDFFile(self.ncfilename,'a')
      ncfile.variables['flags'].assignValue(flags)
      ncfile.history += "\n Flagged on "+time.ctime()
      ncfile.flush()
#      print ncfile.variables['flags'].max()
#      import pdb; pdb.set_trace()
      ncfile.close()
