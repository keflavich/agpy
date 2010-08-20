
import math

import pylab
#from pylab import indices,figure,clf,savefig,plot,legend,text,axes,title,imshow,connect,get_current_fig_manager
from pylab import *
import matplotlib

from mpfit import mpfit

from collapse_gaussfit import *
from ratosexagesimal import *
import pyfits

from numpy import isnan
from mad import MAD,nanmedian

try:
    import coords
except ImportError:
    print "showspec requires coords"

class SpecPlotter:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, cube, axis=None, xtol=None, ytol=None,
          vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x,
          specname=None,dv=None):
    self.vconv = vconv
    self.xtora = xtora
    self.ytodec = ytodec
    self.cube = where(numpy.isnan(cube),0,cube)
    self.specname=specname
    self.dv=dv
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

  def plotspec(self, i, j, fig=None, fignum=1, cube=True,
          button=1, dv=None,ivconv=None,clear=True,color='k',
          **kwargs):
    """
    """
    if dv is None:
        dv = self.dv

    if fig is None and clear:
        fig=figure(fignum)
        fig.clf()
        self.axis = fig.gca()
    elif fig is None:
        self.axis = pylab.gca()
    elif clear:
        fig.clf()
        self.axis = fig.gca()
    #ax = axes([.05,.05,.7,.85])

    vind = self.vconv(arange(self.cube.shape[0]))
    xind = arange(self.cube.shape[0])

    if cube:
        self.axis.plot(vind,self.cube[:,i,j],color=color,
                linestyle='steps-mid',linewidth='.5',
                **kwargs)
    else:
        self.axis.plot(vind,self.cube,color=color,
                linestyle='steps-mid',linewidth='.5',
                **kwargs)
    self.axis.set_xlim(min(vind),max(vind))

    if self.xtora and self.ytodec:
        title("Spectrum at %s %s" % (ratos(self.xtora(i)),dectos(self.ytodec(j))) ) 
    elif self.specname:
        title("Spectrum of %s" % self.specname)
    xlabel("V$_{LSR}$ km s$^{-1}$")
    ylabel("$T_A^*$")
    #legend(loc='best')

def mapplot(plane,cube,vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x):
    
    figure(0)
    clf()
    imshow(plane)

    sp = SpecPlotter(cube,vconv=vconv,xtora=xtora,ytodec=ytodec)
    connect('button_press_event',sp)
 
def splat(filename,vmin=None,vmax=None,button=1,dobaseline=False,exclude=None,
        smooth=None,smoothtype='gaussian',order=1,savepre=None,**kwargs):
    """
    Inputs:
        vmin,vmax - range over which to baseline and plot
        exclude - (internal) range to exclude from baseline fit
    """
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data
    cube = reshape(cube.mean(axis=2).mean(axis=1),[cube.shape[0],1,1])
    dv,v0,p3 = hdr['CD3_3'],hdr['CRVAL3'],hdr['CRPIX3']
    dr,r0,p1 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    dd,d0,p2 = hdr['CD2_2'],hdr['CRVAL2'],hdr['CRPIX2']
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0

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
    if dobaseline: specplot = array([[baseline(cube[argvmin:argvmax].squeeze(),exclude=exclude,order=order)]]).T
    else: specplot = cube[argvmin:argvmax]

    if smooth:
        #specplot[:,0,0] = convolve(specplot[:,0,0],hanning(smooth)/hanning(smooth).sum(),'same')
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot[:,0,0],hanning(smooth)/hanning(smooth).sum(),'same')[::smooth,newaxis,newaxis]
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            specplot = convolve(specplot[:,0,0],kernel,'same')[::smooth,newaxis,newaxis] 
        dv *= smooth
        vconv = lambda v: (v-(p3-argvmin)/smooth+1)*dv+v0
        ivconv = lambda V: (p3-argvmin)/smooth-1+(V-v0)/dv

    sp = SpecPlotter(specplot,vconv=vconv,xtora=xtora,ytodec=ytodec)

    sp.dv = dv
    sp.plotspec(0,0,button=button,ivconv=ivconv,dv=dv,**kwargs)
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

    bspec = baseline(cube,vmin=argvmin,vmax=argvmax,order=order)


def baseline(spectrum,vmin=None,vmax=None,order=1,quiet=True,exclude=None,fitp=None):
    """
    """
    if vmin is None:
        vmin = floor( spectrum.shape[-1]*0.1 )
    if vmax is None:
        vmax = ceil( spectrum.shape[-1]*0.9 )
    
    pguess = [1]*order

    varr = indices(spectrum.shape).squeeze()

    subvarr = varr[vmin:vmax]
    def mpfitfun(data,err):
        def f(p,fjac=None): return [0,numpy.ravel((poly1d(p)(subvarr)-data)/err)]
        return f

    err = ones(spectrum.shape)
    if exclude is not None:
        err[exclude[0]:exclude[1]] = 1e10

    mp = mpfit(mpfitfun(spectrum[vmin:vmax],err[vmin:vmax]),xall=pguess,quiet=quiet)
    fitp = mp.params
    bestfit = poly1d(fitp)(varr).squeeze()

    return (spectrum-bestfit)

def splat_1d(filename,vmin=None,vmax=None,button=1,dobaseline=False,
        exclude=None,smooth=None,order=1,savepre=None,
        smoothtype='hanning',**kwargs):
    """
    """
    f = pyfits.open(filename)
    hdr = f[0].header
    spec = f[0].data
    dv,v0,p3 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    if hdr.get('OBJECT'):
        specname = hdr['OBJECT']
    elif hdr.get('GLON') and hdr.get('GLAT'):
        specname = "%s %s" % (hdr.get('GLON'),hdr.get('GLAT'))
    else:
        specname = filename.remove(".fits")
    if hdr.get('CUNIT1') == 'm/s':
        conversion_factor = 1000.0
    else:
        conversion_factor = 1.0
    vconv = lambda v: ((v-p3+1)*dv+v0)/conversion_factor
    xtora=None
    ytodec=None

    varr = vconv(arange(spec.shape[0]))
    if vmin is None: argvmin = 0
    else: argvmin = argmin(abs(varr-vmin))
    if vmax is None: argvmax = spec.shape[0]
    else: argvmax = argmin(abs(varr-vmax))

    if argvmin > argvmax:
        argvmin,argvmax = argvmax,argvmin
        if exclude is not None: exclude = exclude[::-1]

    if exclude is not None:
        exclude[0] = argmin(abs(varr-exclude[0]))
        exclude[1] = argmin(abs(varr-exclude[1]))
        exclude = array(exclude) - argvmin

    vconv = lambda v: ((v-p3+argvmin+1)*dv+v0) / conversion_factor
    ivconv = lambda V: p3-1-argvmin+(V*conversion_factor-v0)/dv
    if dobaseline: specplot = array([[baseline(spec[argvmin:argvmax].squeeze(),exclude=exclude,order=order)]]).T
    else: specplot = spec[argvmin:argvmax]

    if smooth:
        specplot = convolve(specplot,hanning(smooth)/hanning(smooth).sum(),'same')
    if smooth:
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot,hanning(smooth)/hanning(smooth).sum(),'same')[::smooth]
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            specplot = convolve(specplot,kernel,'same')[::smooth] 
        dv *= smooth
        vconv = lambda v: ((v-(p3-argvmin)/smooth+1)*dv+v0)/conversion_factor
        ivconv = lambda V: (p3-argvmin)/smooth-1+(V*conversion_factor-v0)/dv

    sp = SpecPlotter(specplot,vconv=vconv,xtora=xtora,ytodec=ytodec,specname=specname,dv=dv/conversion_factor)

    sp.plotspec(0,0,button=button,ivconv=ivconv,dv=dv,cube=False,**kwargs)

    if savepre is not None:
        glon,glat = coords.Position([xtora(0),ytodec(0)]).galactic()
        if glat < 0: pm="" 
        else: pm = "+"
        savename = savepre + "G%07.3f%0s%07.3f_" % (glon,pm,glat) + hdr['MOLECULE'].replace(' ','') + hdr['TRANSITI'].replace(' ','')
        savefig(savename+'.png')

    return sp

