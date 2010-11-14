"""
showspec is my homegrown spectrum plotter, meant to somewhat follow STARLINK's
SPLAT and have functionality similar to GAIA, but with an emphasis on producing
publication-quality plots (which, while splat may do, it does unreproducibly)


TO DO:
    -add spectrum arithmetic tools
    -add gaussfitter

"""



import math

import pylab
#from pylab import indices,figure,clf,savefig,plot,legend,text,axes,title,imshow,connect,get_current_fig_manager
from pylab import *
import matplotlib

from mpfit import mpfit

from collapse_gaussfit import *
from ratosexagesimal import *
import pyfits
import gaussfitter

from numpy import isnan
from mad import MAD,nanmedian

try:
    import coords
except ImportError:
    print "showspec requires coords"

def steppify(arr,isX=False,interval=0,sign=+1.0):
    """
    *support function*
    Converts an array to double-length for step plotting
    """
    if isX and interval==0:
        interval = abs(arr[1]-arr[0]) / 2.0
    newarr = array(zip(arr-sign*interval,arr+sign*interval)).ravel()
    return newarr

class SpecPlotter:
    """
    callback for matplotlib to display an annotation when points are clicked on.  The
    point which is closest to the click and within xtol and ytol is identified.
      
    Register this function like this:
      
    scatter(xdata, ydata)
    af = AnnoteFinder(xdata, ydata, annotes)
    connect('button_press_event', af)
    """

    def __init__(self,  cube,  axis=None,  xtol=None,  ytol=None, vconv=lambda x: x, 
            xtora=lambda x: x, ytodec=lambda x: x, specname=None, dv=None,
            hdr=None, errspec=None,  maskspec=None):
      self.vconv = vconv
      self.xtora = xtora
      self.ytodec = ytodec
      self.cube = cube # where(numpy.isnan(cube),0,cube)
      self.spectrum = cube # spectrum is what's plotted; cube is the "raw data"
      self.specname=specname
      self.dv=dv
      self.scale=1.0
      self.units='K'
      self.xunits='km/s'
      self.voff=0.0
      self.offset=0.0
      self.errspec = errspec
      self.erralpha=0.2
      self.plotcolor='k'
      self.specfit = Specfit(self)
      self.baseline = Baseline(self)
      self.vmin=None
      self.vmax=None
      if maskspec is not None:
          self.maskspec = maskspec
      else:
          self.maskspec = zeros(self.cube.shape)
      self.graphic_objects=[]
      if hdr: self.header = hdr
      if axis is None:
        self.axis = pylab.gca()
      else:
        self.axis= axis
  
    def __call__(self, event):
      if event.inaxes:
        clickX = event.xdata
        clickY = event.ydata
        tb = get_current_fig_manager().toolbar
        if ((self.axis is None) or (self.axis==event.inaxes)) and tb.mode=='':
          self.plotspec(clickY,clickX,button=event.button)
  
    def plotspec(self, i=0, j=0, fig=None, fignum=1, cube=False,
            button=1, dv=None,ivconv=None,clear=True,color=None,
            axis=None, offset=None, scale=None, voff=None, vmin=None,
            vmax=None, units=None, xunits=None, erralpha=None, 
            errstyle='fill', **kwargs):
      """
      Plot a spectrum
      Originally written to plot spectra from data cubes, hence the i,j parameter
      to specify the location in the cube
      """
      if dv is None:
          dv = self.dv
  
      if fig is None and clear and axis is None:
          fig=figure(fignum)
          fig.clf()
          self.axis = fig.gca()
      elif fig is None and axis is None:
          self.axis = pylab.gca()
      elif clear and fig:
          fig.clf()
          self.axis = fig.gca()
      elif axis is None:
          self.axis = fig.gca()
      else:
          self.axis = axis
      #ax = axes([.05,.05,.7,.85])
  
      if scale  is not None: self.scale = scale
      if units  is not None: self.units = units
      if xunits is not None: self.xunits= xunits
      if voff   is not None: self.voff  = voff
      if offset is not None: self.offset= offset
      if color  is not None: self.plotcolor=color
      if erralpha is not None: self.erralpha= erralpha
      if vmax  is not None: self.vmax = vmax
      if vmin  is not None: self.vmin = vmin
  
      self.vind = self.vconv(arange(self.cube.shape[0])) + self.voff
      xind = arange(self.cube.shape[0])
  
      if kwargs.has_key('linewidth'):
          linewidth = kwargs.pop('linewidth')
      else:
          linewidth=0.5
  
      if cube:
          self.spectrum = self.cube[:,i,j]*self.scale+self.offset-self.baseline.basespec
          self.axis.plot(self.vind,self.spectrum,color=self.plotcolor,
                  linestyle='steps-mid',linewidth=linewidth,
                  **kwargs)
      else:
          if self.maskspec.sum() > 0:
              nanmask = where(self.maskspec,numpy.nan,1)
              self.spectrum = self.cube*self.scale*nanmask+self.offset-self.baseline.basespec
              self.axis.plot(self.vind,self.spectrum,color=self.plotcolor,
                      linestyle='steps-mid',linewidth=linewidth,
                      **kwargs)
          else:
              self.spectrum = self.cube*self.scale+self.offset-self.baseline.basespec
              self.axis.plot(self.vind,self.spectrum,color=self.plotcolor,
                      linestyle='steps-mid',linewidth=linewidth,
                      **kwargs)
          if self.errspec is not None:
              if errstyle == 'fill':
                  self.axis.fill_between(steppify(self.vind,isX=True,sign=sign(self.dv)),
                          steppify(self.spectrum-self.errspec*self.scale),
                          steppify(self.spectrum+self.errspec*self.scale),
                          facecolor=self.plotcolor, alpha=self.erralpha, **kwargs)
              elif errstyle == 'bars':
                  self.axis.errorbar(self.vind, self.spectrum,
                          yerr=self.errspec*self.scale, ecolor=self.plotcolor, fmt=None,
                          **kwargs)
  
      if vmin is not None: xlo = self.vmin
      else: xlo=self.vind.min()
      if vmax is not None: xhi = self.vmax
      else: xhi=self.vind.max()
      self.axis.set_xlim(xlo,xhi)
  
      if self.xtora and self.ytodec:
          title("Spectrum at %s %s" % (ratos(self.xtora(i)),dectos(self.ytodec(j))) ) 
      elif self.specname:
          title("Spectrum of %s" % self.specname)
      if xunits:
          xlabel(self.xunits)
      else:
          xlabel("V$_{LSR}$ (km s$^{-1}$)")
          self.xunits = 'km/s'
      if units in ['Ta*','Tastar','K']:
        ylabel("$T_A^*$ (K)")
      elif units == 'mJy':
        ylabel("$S_\\nu$ (mJy)")
      elif units == 'Jy':
        ylabel("$S_\\nu$ (Jy)")
      else:
        ylabel(self.units)
      #legend(loc='best')
  
    def save(self,fname,**kwargs):
      """
      Save the current spectrum (useful for saving baselined data)
      """
      newfile = pyfits.PrimaryHDU(data=self.cube,header=self.header)
      newfile.writeto(fname,**kwargs)
  
    def showlines(self,linefreqs,linenames,ctype='freq',cunit='hz',yscale=0.8,vofflines=0.0,
            voffunit='km/s',**kwargs):
        """
        Overplot vertical lines and labels at the frequencies (or velocities) of each line
  
        yscale - fraction of maximum at which to label
        """
  
        if ctype != 'freq':
            print "Sorry, non-frequency units not implemented yet."
            return
  
        speedoflight=2.99792458e5
        if 'hz' in cunit or 'Hz' in cunit:
            linefreqs *= (1.0 + vofflines / speedoflight)
        else:
            linefreqs += vofflines
      
        ymax = (self.spectrum[self.spectrum==self.spectrum]).max() 
        for lf,ln in zip(lineh2cofreqs,linenames):
            if lf < self.vind.max() and lf > self.vind.min():
                self.graphic_objects.append(vlines(lf,0,ymax,**kwargs))
                self.graphic_objects.append(text(lf,ymax*yscale,ln,rotation='vertical',**kwargs))
  
        return self.graphic_objects
  
    def hidelines(self):
          for obj in self.graphic_objects:
              obj.set_visible(False)
  
    def fitspec(self, interactive=False, usemoments=True, fitcolor='r',
            multifit=False, guesses=None, annotate=True, 
            **kwargs):
        """
        Fit gaussians to a spectrum
        """
  
        self.specfit.fitcolor = fitcolor
        self.specfit.clear()

        self.ngauss = 0
        self.specfit.fitkwargs = kwargs
        if interactive:
            print "Left-click twice to select a fitting range, then middle-click twice to select a peak and width"
            self.specfit.nclicks_b1 = 0
            self.specfit.nclicks_b2 = 0
            self.specfit.guesses = []
            self.specfit.click = self.axis.figure.canvas.mpl_connect('button_press_event',self.makeguess)
            self.specfit.autoannotate = annotate
        elif multifit:
            if guesses is None:
                print "You must input guesses when using multifit.  Also, baseline first!"
            else:
                self.specfit.guesses = guesses
                self.specfit.multifit()
                self.specfit.autoannotate = annotate
        else:
            print "Non-interactive, 1D fit with automatic guessing"
            if self.baseline.order is None:
                self.specfit.onedfit(usemoments=usemoments,annotate=annotate)
            else:
                self.specfit.onedfit(usemoments=usemoments,annotate=annotate,
                        vheight=False,height=0.0)
            self.refresh()
    
    def makeguess(self,event):
        if event.button == 1:
            self.specfit.selectregion(event)
        elif event.button == 2:
            self.specfit.guesspeakwidth(event)
        elif event.button == 3:
            disconnect(self.specfit.click)
            if self.specfit.ngauss > 0:
                print len(self.specfit.guesses)/3," Guesses: ",self.specfit.guesses," X channel range: ",self.specfit.gx1,self.specfit.gx2
                if len(self.specfit.guesses) % 3 == 0:
                    self.specfit.multifit()
                    for p in self.specfit.guessplot + self.specfit.fitregion:
                        p.set_visible(False)
                else: 
                    print "error, wrong # of pars"
        self.refresh()

    def refresh(self):
        self.axis.figure.canvas.draw()

class Baseline:
    def __init__(self,specplotter):
        self.baselinepars  = None
        self.order = None
        self.basespec = zeros(specplotter.spectrum.shape[0])
        self.specplotter = specplotter
        self.blleg = None

    def __call__(self,order=1,annotate=False,**kwargs):
        fitp = zeros(order+1)
        self.specplotter.spectrum, self.baselinepars = baseline(
                self.specplotter.spectrum,
                xarr=self.specplotter.vind,
                order=order, **kwargs)
        self.order = order
        self.basespec = poly1d(self.baselinepars)(self.specplotter.vind)
        self.specplotter.plotspec()
        if annotate: self.annotate()

    def annotate(self,xloc=0.02,yloc=0.95):
        bltext = "bl:"+"".join(["$%gx^{%i}$" % (f,i) for i,f in enumerate(self.baselinepars)])
        #self.blleg = text(xloc,yloc     ,bltext,transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        self.blleg = legend((self.specplotter.axis.lines[0],),(bltext,),loc='upper left',markerscale=0)
        self.specplotter.axis.add_artist(self.blleg)
        self.specplotter.refresh()
  
    def clearlegend(self):
        if self.blleg is not None: 
            self.blleg.set_visible(False)
            self.specplotter.axis.artists.remove(self.blleg)

class Specfit:

    def __init__(self,specplotter):
        self.model = None
        self.modelpars = None
        self.modelplot = None
        self.guessplot = []
        self.fitregion = []
        self.ngauss = 0
        self.nclicks_b1 = 0
        self.nclicks_b2 = 0
        self.gx1 = 0
        self.gx2 = specplotter.spectrum.shape[0]
        self.guesses = []
        self.click = 0
        self.fitkwargs = {}
        self.auto = 0
        self.autoannotate = True
        self.specplotter = specplotter
        self.gaussleg=None

        if self.specplotter.errspec is not None:
            self.errspec = self.specplotter.errspec[self.gx1:self.gx2]
        else: self.errspec = ones(specplotter.spectrum.shape[0])

    def multifit(self):
  
        mpp,model,mpperr,chi2 = gaussfitter.multigaussfit(
                self.specplotter.vind[self.gx1:self.gx2], 
                self.specplotter.spectrum[self.gx1:self.gx2], 
                err=self.errspec,
                ngauss=len(self.guesses)/3,
                params=self.guesses,
                **self.fitkwargs)
        self.model = model
        self.modelpars = mpp
        self.modelplot = self.specplotter.axis.plot(
                self.specplotter.vind[self.gx1:self.gx2],
                self.model, color=self.fitcolor, linewidth=0.5)
        if self.autoannotate:
            self.annotate()
    
    def onedfit(self, usemoments=True, annotate=True, vheight=True, height=0):
        self.ngauss += 1
        self.auto = 1
        mpp,model,mpperr,chi2 = gaussfitter.onedgaussfit(
                self.specplotter.vind[self.gx1:self.gx2],
                self.specplotter.spectrum[self.gx1:self.gx2],
                err=self.errspec,
                vheight=vheight,
                params=[height,1,0,1],
                usemoments=usemoments,
                **self.fitkwargs)
        self.model = model
        self.modelpars = mpp[1:]
        self.specplotter.baseline.baselinepars = mpp[:1]
        self.modelplot = self.specplotter.axis.plot(
                self.specplotter.vind[self.gx1:self.gx2],
                self.model, color=self.fitcolor, linewidth=0.5)
        if annotate:
            self.annotate()
            if vheight: self.specplotter.baseline.annotate()

    def plotresiduals(self,figure=figure(2)):
        self.residualfigure = figure
        self.residualplot = self.residualfigure.axis.plot(self.specplotter.vind[self.gx1:self.gx2],
                self.specplotter.cube[self.gx1:self.gx2]-self.model,drawstyle='steps-mid',
                linewidth=0.5, color='k')

    def annotate(self,xloc=0.8,yloc=0.95):
        #text(xloc,yloc     ,"c=%g" % self.modelpars[1],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.05,"w=%g" % self.modelpars[2],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.10,"a=%g" % self.modelpars[0],transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        self.gaussleg = legend((self.specplotter.axis.lines[0],self.specplotter.axis.lines[0],self.specplotter.axis.lines[0]),
                ("c=%g" % self.modelpars[1],
                "w=%g" % self.modelpars[2],
                "a=%g" % self.modelpars[0]),
                loc='upper right',markerscale=0)
        self.specplotter.axis.add_artist(self.gaussleg)
        self.specplotter.refresh()

    def selectregion(self,event):
        if self.nclicks_b1 == 0:
            self.gx1 = argmin(abs(event.xdata-self.specplotter.vind))
            self.nclicks_b1 += 1
        elif self.nclicks_b1 == 1:
            self.gx2 = argmin(abs(event.xdata-self.specplotter.vind))
            self.nclicks_b1 -= 1
            if self.gx1 > self.gx2: self.gx1,self.gx2 = self.gx2,self.gx1
            self.fitregion = self.specplotter.axis.plot(
                    self.specplotter.vind[self.gx1:self.gx2],
                    self.specplotter.spectrum[self.gx1:self.gx2],
                    drawstyle='steps-mid',
                    color='c')
            if self.guesses == []:
                self.guesses = gaussfitter.onedmoments(
                        self.specplotter.vind[self.gx1:self.gx2],
                        self.specplotter.spectrum[self.gx1:self.gx2],
                        vheight=0)
                self.ngauss = 1
                self.auto = 1
            if self.specplotter.errspec is not None:
                self.errspec = self.specplotter.errspec[self.gx1:self.gx2]
            else: self.errspec = ones(self.gx2-self.gx1)

    def guesspeakwidth(self,event):
        if self.nclicks_b2 == 0:
            if (self.auto == 1 and self.ngauss == 1):
                self.guesses[:2] = [event.ydata,event.xdata]
            else:
                self.guesses += [event.ydata,event.xdata,1]
                self.ngauss += 1
            self.nclicks_b2 += 1
            self.guessplot = [self.specplotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
        elif self.nclicks_b2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2])
            self.nclicks_b2 -= 1
            self.guessplot += self.specplotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')

    def clear(self):
        if self.modelplot is not None:
            for p in self.modelplot:
                p.set_visible(False)

    def clearlegend(self):
        if self.gaussleg is not None: 
            self.gaussleg.set_visible(False)
            self.specplotter.axis.artists.remove(self.gaussleg)
        

def mapplot(plane,cube,vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x):
    
    figure(0)
    clf()
    imshow(plane)

    sp = SpecPlotter(cube,vconv=vconv,xtora=xtora,ytodec=ytodec)
    connect('button_press_event',sp)


def open_3d(filename):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data
    if len(cube.shape) == 4: cube=cube[0,:,:,:]
    cube = reshape(cube.mean(axis=2).mean(axis=1),[cube.shape[0],1,1])
    dv,v0,p3 = hdr['CD3_3'],hdr['CRVAL3'],hdr['CRPIX3']
    dr,r0,p1 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    dd,d0,p2 = hdr['CD2_2'],hdr['CRVAL2'],hdr['CRPIX2']
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0

    return dv,v0,p3,hdr,cube,xtora,ytodec,vconv
 
def splat(filename,vmin=None,vmax=None,button=1,dobaseline=False,exclude=None,
        smooth=None,smoothtype='gaussian',order=1,savepre=None,**kwargs):
    """
    Inputs:
        vmin,vmax - range over which to baseline and plottransform = ax.transAxes
        exclude - (internal) range to exclude from baseline fit
    """
    dv,v0,p3,hdr,cube,xtora,ytodec,vconv = open_3d(filename)

    splat_1d(vpars=[dv,v0,p3],hdr=hdr,cube=cube[:,0,0],xtora=xtora,ytodec=ytodec,vconv=vconv)

    varr = vconv(arange(cube.shape[0]))
    if vmin is None: argvmin = 0
    else: argvmin = argmin(abs(varr-vmin))
    if vmax is None: argvmax = cube.shape[0]
    else: argvmax = argmin(abs(varr-vmax))

    if argvmin > argvmax:
        argvmin,argvmax = argvmax,argvmin
        if exclude is not None: exclude = exclude[::-1]

    if exclude is not None:
        exclude[0] = argmin(abs(varr-exclude[0]))
        exclude[1] = argmin(abs(varr-exclude[1]))
        exclude = array(exclude) - argvmin

    vconv = lambda v: (v-p3+argvmin+1)*dv+v0
    ivconv = lambda V: p3-1-argvmin+(V-v0)/dv
    #if dobaseline: specplot = array([[baseline(cube[argvmin:argvmax].squeeze(),exclude=exclude,order=order)]]).T
    specplot = cube[argvmin:argvmax]

    if smooth:
        #specplot[:,0,0] = convolve(specplot[:,0,0],hanning(smooth)/hanning(smooth).sum(),'same')
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot[:,0,0],hanning(2+smooth)/hanning(2+smooth).sum(),'same')[::smooth,newaxis,newaxis]
        elif smoothtype == 'boxcar':
            specplot = convolve(specplot[:,0,0],ones(smooth)/float(smooth),'same')[::smooth,newaxis,newaxis]
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            specplot = convolve(specplot[:,0,0],kernel,'same')[::smooth,newaxis,newaxis] 
        # this bit of code may also make sense, but I'm shifting the center pixel instead
        # b/c it's easier (?) to deal with velocity range
        #v0 += (abs(dv)*smooth - abs(dv))/2.0 # pixel center moves by half the original pixel size
        dv *= smooth
        newrefpix = (p3-0.5-argvmin)/smooth  # this was resolved by advanced guess-and check
        # but also, sort of makes sense: FITS refers to the *center* of a pixel.  You want to 
        # shift 1/2 pixel to the right so that the first pixel goes from 0 to 1
        vconv = lambda v: ((v-newrefpix)*dv+v0)/conversion_factor
        ivconv = lambda V: newrefpix+(V*conversion_factor-v0)/dv

    sp = SpecPlotter(specplot,vconv=vconv,xtora=xtora,ytodec=ytodec)

    sp.dv = dv
    sp.plotspec(button=button,ivconv=ivconv,dv=dv,cube=True,**kwargs)
    if dobaseline: 
        sp.baseline(exclude=exclude,order=order,quiet=quiet)
    sp.axis.set_xlim(vmin,vmax)

    if savepre is not None:
        glon,glat = coords.Position([xtora(0),ytodec(0)]).galactic()
        if glat < 0: pm="" 
        else: pm = "+"
        savename = savepre + "G%07.3f%0s%07.3f_" % (glon,pm,glat) + hdr['MOLECULE'].replace(' ','') + hdr['TRANSITI'].replace(' ','')
        savefig(savename+'.png')

    return sp

def gaia(filename,estimator='max',axis=0):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data
    dv,v0,p3 = hdr['CD3_3'],hdr['CRVAL3'],hdr['CRPIX3']
    dr,r0,p1 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    dd,d0,p2 = hdr['CD2_2'],hdr['CRVAL2'],hdr['CRPIX2']
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0

    if axis > 0:
        cube = cube.swapaxes(0,axis)

    if estimator == 'max':
        p = where(isnan(cube),0,cube).max(axis=0)
    elif estimator == 'int':
        p = where(isnan(cube),0,cube).sum(axis=0) * dv
    elif estimator == 'intdivmax':
        cut = MAD(cube.ravel()) + nanmedian(cube.ravel())
        if cut < 0:
            cut = 0
        m = where(isnan(cube),0,cube).max(axis=0)
        i = where(isnan(cube),0,cube).sum(axis=0) * dv
        p = where(i<0,0,i)/where(m<=cut,numpy.inf,m)
    elif estimator[-5:] == ".fits":
        p = pyfits.open(estimator)[0].data

    mapplot(p,cube,vconv,xtora,ytodec)

def baseline_file(filename,outfilename,vmin=None,vmax=None,order=1,crop=False):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data.squeeze()
    dv,v0,p3 = hdr['CD3_3'],hdr['CRVAL3'],hdr['CRPIX3']
    dr,r0,p1 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    dd,d0,p2 = hdr['CD2_2'],hdr['CRVAL2'],hdr['CRPIX2']
    vconv = lambda v: (v-p3+1)*dv+v0
    varr = vconv(arange(cube.shape[-1]))
    if vmin is None: argvmin = None
    else: argvmin = argmin(abs(varr-vmin))
    if vmax is None: argvmax = None
    else: argvmax = argmin(abs(varr-vmax))

    bspec,bfit = baseline(cube,vmin=argvmin,vmax=argvmax,order=order)


def baseline(spectrum,xarr=None,xmin=None,xmax=None,order=1,quiet=True,exclude=None):
    """
    Subtract a baseline from a spectrum
    If xmin,xmax are not specified, defaults to ignoring first and last 10% of spectrum

    exclude is a set of start/end indices to ignore when baseline fitting
    (ignored by setting error to infinite in fitting procedure)
    """
    if xmin is None:
        xmin = floor( spectrum.shape[-1]*0.1 )
    if xmax is None:
        xmax = ceil( spectrum.shape[-1]*0.9 )
    
    pguess = [1]*(order+1)

    if xarr is None:
        xarr = indices(spectrum.shape).squeeze()

    subxarr = xarr[xmin:xmax]
    def mpfitfun(data,err):
        def f(p,fjac=None): return [0,numpy.ravel((poly1d(p)(subxarr)-data)/err)]
        return f

    err = ones(spectrum.shape)
    if exclude is not None:
        err[exclude[0]:exclude[1]] = 1e10

    mp = mpfit(mpfitfun(spectrum[xmin:xmax],err[xmin:xmax]),xall=pguess,quiet=quiet)
    fitp = mp.params
    bestfit = poly1d(fitp)(xarr).squeeze()

    return (spectrum-bestfit),fitp

def open_1d(filename,specnum=0,wcstype='',errspecnum=None,maskspecnum=None):
    """
    Grabs all the relevant pieces of a 1d spectrum for plotting
    wcstype is the suffix on the WCS type to get to velocity/frequency/whatever
    """
    f = pyfits.open(filename)
    hdr = f[0].header
    spec = f[0].data
    errspec  = None
    maskspec = None
    if hdr.get('NAXIS') == 2:
        if errspecnum is not None:
            errspec = spec[errspecnum,:]
        if maskspecnum is not None:
            maskspec = spec[maskspecnum,:]
        spec = spec[specnum,:]
    elif hdr.get('NAXIS') > 2:
        raise ValueError("Too many axes for open_1d (splat_1d) - use cube instead")
    if hdr.get('CD1_1'+wcstype):
        dv,v0,p3 = hdr['CD1_1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    else:
        dv,v0,p3 = hdr['CDELT1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    if hdr.get('OBJECT'):
        specname = hdr['OBJECT']
    elif hdr.get('GLON') and hdr.get('GLAT'):
        specname = "%s %s" % (hdr.get('GLON'),hdr.get('GLAT'))
    else:
        specname = filename.rstrip(".fits")
    if hdr.get('CUNIT1'+wcstype) in ['m/s','M/S']:
        conversion_factor = 1000.0
        xunits = 'm/s'
    else:
        conversion_factor = 1.0
        xunits = hdr.get('CUNIT1'+wcstype)
    vconv = lambda v: ((v-p3+1)*dv+v0)/conversion_factor
    xtora=None
    ytodec=None
    units = hdr.get('BUNIT')
    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'

    return dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname,units,xunits,errspec,maskspec

def splat_1d(filename=None,vmin=None,vmax=None,button=1,dobaseline=False,
        exclude=None,smooth=None,order=1,savepre=None,vcrop=True,
        vconv=None,vpars=None,hdr=None,spec=None,xtora=None,ytora=None,
        specname=None,quiet=True,specnum=0,errspecnum=None,wcstype='',
        offset=None, annotatebaseline=True,
        smoothtype='gaussian',convmode='valid',maskspecnum=None,**kwargs):
    """
    Inputs:
        vmin,vmax - range over which to baseline and plot
        exclude - (internal) range to exclude from baseline fit
        vcrop - will vmin/vmax crop out data, or just set the plot limits?
    """
    if vpars and vconv and hdr and spec and xtora and ytora:
        dv,v0,p3 = vpars
    else:
        dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname_file,units,xunits,errspec,maskspec = \
                open_1d(filename,specnum=specnum,wcstype=wcstype,errspecnum=errspecnum,maskspecnum=maskspecnum)
        if specname is None: specname=specname_file
        if units is None and kwargs.has_key('units'): units = kwargs.pop('units')

    if type(offset)==type('str'):
        if hdr.get(offset) is not None:
            offset = hdr.get(offset)
        else:
            raise ValueError("Offset specified but none present.")
    elif offset is None:
        offset = 0.0

    varr = vconv(arange(spec.shape[0]))
    if vmin is None or vcrop==False: argvmin = 0
    else: 
      argvmin = argmin(abs(varr-vmin))
      if dv > 0:
        hdr.update('CRPIX1'+wcstype,p3-argvmin)
    if vmax is None or vcrop==False: argvmax = spec.shape[0]
    else: 
      argvmax = argmin(abs(varr-vmax))
      if dv < 0:
        hdr.update('CRPIX1'+wcstype,p3-argvmax)

    if argvmin > argvmax:
        argvmin,argvmax = argvmax,argvmin
        if exclude is not None: exclude = exclude[::-1]

    if exclude is not None:
        exclude[0] = argmin(abs(varr-exclude[0]))
        exclude[1] = argmin(abs(varr-exclude[1]))
        exclude = array(exclude) - argvmin

    vconv = lambda v: ((v-p3+argvmin+1)*dv+v0) / conversion_factor
    ivconv = lambda V: p3-1-argvmin+(V*conversion_factor-v0)/dv

    specplot = spec[argvmin:argvmax]
    if errspec is not None: errspec=errspec[argvmin:argvmax]
    if maskspec is not None: maskspec=maskspec[argvmin:argvmax]

    if smooth:
        roundsmooth = round(smooth) # can only downsample by integers
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot,hanning(2+roundsmooth)/hanning(2+roundsmooth).sum(),convmode)[::roundsmooth]
            kernsize = smooth
            ones_sameshape = zeros(smooth+2)
            ones_sameshape[1:-1] = 1
        elif smoothtype == 'boxcar':
            specplot = convolve(specplot,ones(roundsmooth)/float(roundsmooth),convmode)[::roundsmooth]
            kernsize = roundsmooth
            ones_sameshape = ones(roundsmooth)
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            kernsize = len(kernel)
            specplot = convolve(specplot,kernel,convmode)[::roundsmooth] 
            ones_sameshape = zeros(roundsmooth*3)
            ones_sameshape[roundsmooth:-roundsmooth] = 1
        if errspec is not None: 
            errspec = sqrt(convolve(errspec**2,ones_sameshape,convmode)[::roundsmooth]) / float(roundsmooth)
        if maskspec is not None: 
            maskspec = array(convolve(maskspec,ones_sameshape,convmode)[::roundsmooth],dtype='bool')
            if maskspec.shape != specplot.shape: import pdb; pdb.set_trace()
        # this bit of code may also make sense, but I'm shifting the center pixel instead
        # b/c it's easier (?) to deal with velocity range
        #v0 += (abs(dv)*smooth - abs(dv))/2.0 # pixel center moves by half the original pixel size
        dv *= roundsmooth
        if convmode == 'same':
          newrefpix = (p3-argvmin)/roundsmooth  
        elif convmode == 'full':
          newrefpix = (p3-0.5-argvmin+kernsize/2.0)/roundsmooth  
        elif convmode == 'valid':
          newrefpix = (p3-0.5-argvmin-kernsize/2.0)/roundsmooth  
        # this was resolved by advanced guess-and check
        # but also, sort of makes sense: FITS refers to the *center* of a pixel.  You want to 
        # shift 1/2 pixel to the right so that the first pixel goes from 0 to 1
        vconv = lambda v: ((v-newrefpix)*dv+v0)/conversion_factor
        ivconv = lambda V: newrefpix+(V*conversion_factor-v0)/dv
        hdr.update('CRPIX1'+wcstype,newrefpix+1)
        hdr.update('CDELT1'+wcstype,dv)

    sp = SpecPlotter(specplot, vconv=vconv, xtora=xtora, ytodec=ytodec,
            specname=specname, dv=dv/conversion_factor, hdr=hdr,
            errspec=errspec, maskspec=maskspec)

    sp.plotspec(button=button, ivconv=ivconv, dv=dv, cube=False,
            vmin=vmin, vmax=vmax, units=units, xunits=xunits, offset=offset,
            **kwargs)
    if dobaseline: 
        sp.baseline(exclude=exclude,order=order,quiet=quiet,annotate=annotatebaseline)
    sp.refresh()
    
    if hdr.get('GLON') and hdr.get('GLAT'):
        sp.glon = hdr.get('GLON')
        sp.glat = hdr.get('GLAT')

    if savepre is not None:
        glon,glat = sp.glon,sp.glat
        if glat < 0: pm="" 
        else: pm = "+"
        savename = savepre + "G%07.3f%0s%07.3f_" % (glon,pm,glat) + hdr['MOLECULE'].replace(' ','') + hdr['TRANSITI'].replace(' ','')
        savefig(savename+'.png')

    return sp

def splat_tspec(filename,specnum=0,**kwargs):
    """
    Same as splat_1d for tspec data
    """

    tdata = pyfits.getdata(filename)
    theader = pyfits.getheader(filename)
    if len(tdata.shape) == 3:
        tdata = tdata[specnum,:,:]
    wavelength = tdata[0,:]
    spectrum   = tdata[1,:]
    error      = tdata[2,:]

    vconv = lambda x: wavelength[x]
    ivconv = lambda x: argmin(abs(wavelength-x))
    
    specname='TSPEC'
    dv = median(wavelength[1:] - wavelength[:-1])
    
    sp = SpecPlotter(spectrum,vconv=vconv,specname=specname,dv=dv,hdr=theader)

    sp.plotspec(ivconv=ivconv,dv=dv,cube=False,units=theader.get('YUNITS'),xunits=theader.get('XUNITS'),**kwargs)

    return sp
