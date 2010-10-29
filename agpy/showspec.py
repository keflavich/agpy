
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
          specname=None,dv=None,hdr=None):
    self.vconv = vconv
    self.xtora = xtora
    self.ytodec = ytodec
    self.cube = where(numpy.isnan(cube),0,cube)
    self.specname=specname
    self.dv=dv
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

  def plotspec(self, i, j, fig=None, fignum=1, cube=True,
          button=1, dv=None,ivconv=None,clear=True,color='k',
          axis=None, offset=0.0, scale=1.0, voff=0.0, vmin=None,
          vmax=None, units='K', xunits=None, **kwargs):
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

    vind = self.vconv(arange(self.cube.shape[0])) + voff
    xind = arange(self.cube.shape[0])

    if kwargs.has_key('linewidth'):
        linewidth = kwargs.pop('linewidth')
    else:
        linewidth="0.5"

    if cube:
        self.axis.plot(vind,self.cube[:,i,j]*scale+offset,color=color,
                linestyle='steps-mid',linewidth=linewidth,
                **kwargs)
    else:
        self.axis.plot(vind,self.cube*scale+offset,color=color,
                linestyle='steps-mid',linewidth=linewidth,
                **kwargs)

    if vmin is not None: xlo = vmin
    else: xlo=vind.min()
    if vmax is not None: xhi = vmax
    else: xhi=vind.max()
    self.axis.set_xlim(xlo,xhi)

    if self.xtora and self.ytodec:
        title("Spectrum at %s %s" % (ratos(self.xtora(i)),dectos(self.ytodec(j))) ) 
    elif self.specname:
        title("Spectrum of %s" % self.specname)
    if xunits:
        xlabel(xunits)
    else:
        xlabel("V$_{LSR}$ (km s$^{-1}$)")
    self.units = units
    if units in ['Ta*','Tastar','K']:
      ylabel("$T_A^*$ (K)")
    elif units == 'mJy':
      ylabel("$S_\\nu$ (mJy)")
    elif units == 'Jy':
      ylabel("$S_\\nu$ (Jy)")
    else:
      ylabel(units)
    #legend(loc='best')

  def save(self,fname,**kwargs):
    """
    Save the current spectrum (useful for saving baselined data)
    """
    newfile = pyfits.PrimaryHDU(data=self.cube,header=self.header)
    newfile.writeto(fname,**kwargs)

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
        vmin,vmax - range over which to baseline and plot
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
    if dobaseline: specplot = array([[baseline(cube[argvmin:argvmax].squeeze(),exclude=exclude,order=order)]]).T
    else: specplot = cube[argvmin:argvmax]

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
    Subtract a baseline from a spectrum
    If vmin,vmax are not specified, defaults to ignoring first and last 10% of spectrum

    exclude is a set of start/end indices to ignore when baseline fitting
    (ignored by setting error to infinite in fitting procedure)
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

def open_1d(filename,specnum=0):
    f = pyfits.open(filename)
    hdr = f[0].header
    spec = f[0].data
    if hdr.get('NAXIS') == 2:
        spec = spec[specnum,:]
    elif hdr.get('NAXIS') > 2:
        raise ValueError("Too many axes for open_1d (splat_1d) - use cube instead")
    if hdr.get('CD1_1'):
        dv,v0,p3 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
    else:
        dv,v0,p3 = hdr['CDELT1'],hdr['CRVAL1'],hdr['CRPIX1']
    if hdr.get('OBJECT'):
        specname = hdr['OBJECT']
    elif hdr.get('GLON') and hdr.get('GLAT'):
        specname = "%s %s" % (hdr.get('GLON'),hdr.get('GLAT'))
    else:
        specname = filename.rstrip(".fits")
    if hdr.get('CUNIT1') in ['m/s','M/S']:
        conversion_factor = 1000.0
    else:
        conversion_factor = 1.0
    vconv = lambda v: ((v-p3+1)*dv+v0)/conversion_factor
    xtora=None
    ytodec=None
    units = hdr.get('BUNIT')

    return dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname,units

def splat_1d(filename=None,vmin=None,vmax=None,button=1,dobaseline=False,
        exclude=None,smooth=None,order=1,savepre=None,vcrop=True,
        vconv=None,vpars=None,hdr=None,spec=None,xtora=None,ytora=None,
        specname=None,quiet=True,specnum=0,
        smoothtype='gaussian',convmode='valid',**kwargs):
    """
    Inputs:
        vmin,vmax - range over which to baseline and plot
        exclude - (internal) range to exclude from baseline fit
        vcrop - will vmin/vmax crop out data, or just set the plot limits?
    """
    if vpars and vconv and hdr and spec and xtora and ytora:
        dv,v0,p3 = vpars
    else:
        dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname_file,units = open_1d(filename,specnum=specnum)
        if specname is None: specname=specname_file

    varr = vconv(arange(spec.shape[0]))
    if vmin is None or vcrop==False: argvmin = 0
    else: 
      argvmin = argmin(abs(varr-vmin))
      if dv > 0:
        hdr.update('CRPIX1',p3-argvmin)
    if vmax is None or vcrop==False: argvmax = spec.shape[0]
    else: 
      argvmax = argmin(abs(varr-vmax))
      if dv < 0:
        hdr.update('CRPIX1',p3-argvmax)

    if argvmin > argvmax:
        argvmin,argvmax = argvmax,argvmin
        if exclude is not None: exclude = exclude[::-1]

    if exclude is not None:
        exclude[0] = argmin(abs(varr-exclude[0]))
        exclude[1] = argmin(abs(varr-exclude[1]))
        exclude = array(exclude) - argvmin

    vconv = lambda v: ((v-p3+argvmin+1)*dv+v0) / conversion_factor
    ivconv = lambda V: p3-1-argvmin+(V*conversion_factor-v0)/dv
    if dobaseline: specplot = baseline(spec[argvmin:argvmax].squeeze(),exclude=exclude,order=order,quiet=quiet)
    else: specplot = spec[argvmin:argvmax]

    if smooth:
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot,hanning(2+smooth)/hanning(2+smooth).sum(),convmode)[::smooth]
            kernsize = smooth
        elif smoothtype == 'boxcar':
            specplot = convolve(specplot,ones(smooth)/float(smooth),convmode)[::smooth]
            kernsize = smooth
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            kernsize = len(kernel)
            specplot = convolve(specplot,kernel,convmode)[::smooth] 
        # this bit of code may also make sense, but I'm shifting the center pixel instead
        # b/c it's easier (?) to deal with velocity range
        #v0 += (abs(dv)*smooth - abs(dv))/2.0 # pixel center moves by half the original pixel size
        dv *= smooth
        if convmode == 'same':
          newrefpix = (p3-argvmin)/smooth  
        elif convmode == 'full':
          newrefpix = (p3-0.5-argvmin+kernsize/2.0)/smooth  
        elif convmode == 'valid':
          newrefpix = (p3-0.5-argvmin-kernsize/2.0)/smooth  
        # this was resolved by advanced guess-and check
        # but also, sort of makes sense: FITS refers to the *center* of a pixel.  You want to 
        # shift 1/2 pixel to the right so that the first pixel goes from 0 to 1
        vconv = lambda v: ((v-newrefpix)*dv+v0)/conversion_factor
        ivconv = lambda V: newrefpix+(V*conversion_factor-v0)/dv
        hdr.update('CRPIX1',newrefpix+1)
        hdr.update('CDELT1',dv)

    sp = SpecPlotter(specplot,vconv=vconv,xtora=xtora,ytodec=ytodec,specname=specname,dv=dv/conversion_factor,hdr=hdr)

    sp.plotspec(0,0,button=button,ivconv=ivconv,dv=dv,cube=False,vmin=vmin,vmax=vmax,units=units,**kwargs)
    
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

    sp.plotspec(0,0,ivconv=ivconv,dv=dv,cube=False,units=theader.get('YUNITS'),xunits=theader.get('XUNITS'),**kwargs)

    return sp
