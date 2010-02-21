
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

  def __init__(self, cube, axis=None, xtol=None, ytol=None,vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x):
    self.vconv = vconv
    self.xtora = xtora
    self.ytodec = ytodec
    self.cube = where(numpy.isnan(cube),0,cube)
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

  def plotspec(self, i, j, fignum=1, button=1, dv=None,ivconv=None,
          vguess1=None,vguess2=None,vguess3=None,wguess1=None,wguess2=None,wguess3=None,
          screenprint=True,fprint=None):
    """
    """
    f=figure(fignum)
    f.clf()
    ax = axes([.05,.05,.7,.85])

    vind = self.vconv(arange(self.cube[:,0,0].shape[0]))
    xind = arange(self.cube[:,0,0].shape[0])

    pars = return_param(self.cube[:,i,j])
    doublepars = return_double_param(self.cube[:,i,j])
    if vguess1 is None: vguess1 = (doublepars[0])
    elif ivconv is not None: vguess1=ivconv(vguess1)
    if vguess2 is None: vguess2 = (doublepars[1])
    elif ivconv is not None: vguess2=ivconv(vguess2)
    if vguess3 is None: vguess3 = ((doublepars[0]+doublepars[1])/2.)
    elif ivconv is not None: vguess3=ivconv(vguess3)
    if wguess1 is None: wguess1 = doublepars[2]    /dv
    if wguess2 is None: wguess2 = doublepars[3]    /dv
    if wguess3 is None: wguess3 = doublepars[2]*5. /dv
    tpguess = [vguess1,vguess2,vguess3,wguess1,wguess2,wguess3,doublepars[4],doublepars[5],doublepars[4]/5.]
    triplepars = return_triple_param(self.cube[:,i,j],params=tpguess)

    plot(vind,self.cube[:,i,j],color='black',linestyle='steps',linewidth='.5')

    chi2_1 = (((gerr(self.cube[:,i,j])(pars))**2).sum() )
    chi2_2 = ((double_gerr(self.cube[:,i,j])(doublepars))**2).sum()
    chi2_3 = ((triple_gerr(self.cube[:,i,j])(triplepars))**2).sum()

    if button==1:
        plot(vind,gaussian(*pars)(xind),'r-.',label="Single %f" % ( chi2_1 ))
        plot(vind,double_gaussian(*doublepars)(xind),'g:',label="Double %f" % ( chi2_2 ))
    elif button==2:
        plot(vind,gaussian(triplepars[0],triplepars[3],triplepars[6])(xind),'r:')
        plot(vind,gaussian(triplepars[1],triplepars[4],triplepars[7])(xind),'r:')
        plot(vind,gaussian(triplepars[2],triplepars[5],triplepars[8])(xind),'r:')

    plot(vind,triple_gaussian(*triplepars)(xind),'b--',label="Triple %f" % ( chi2_3 ),linewidth=2)

    pars[0] = self.vconv(pars[0])
    text(1.05,.8,"c1 %3.2f w1 %3.2f a1 %3.2f" % tuple(pars),transform=ax.transAxes,size='smaller')
    dp = [ self.vconv(doublepars[0]) , doublepars[2]/dv, doublepars[4],
            self.vconv(doublepars[1]), doublepars[3]/dv, doublepars[5] ]
    text(1.05,.6,"c1 %3.2f w1 %3.2f a1 %3.2f\nc2 %3.2f w2 %3.2f a2 %3.2f" % tuple(dp),transform=ax.transAxes,size='smaller')
    tp = [ self.vconv(triplepars[0]) , triplepars[3]/dv, triplepars[6],
            self.vconv(triplepars[1]), triplepars[4]/dv, triplepars[7],
            self.vconv(triplepars[2]), triplepars[5]/dv, triplepars[8]  ]
    text(1.05,.4,"c1 %3.2f w1 %3.2f a1 %3.2f\nc2 %3.2f w2 %3.2f a2 %3.2f\nc3 %3.2f w3 %3.2f a3 %3.2f" % tuple(tp),transform=ax.transAxes,size='smaller')

    if screenprint:
        print "%12.2f%12.2f%3.2f" % tuple(pars) + \
                "%12.2f%12.2f%12.2f  %12.2f%12.2f%3.2f" % tuple(dp) + \
                "%12.2f%12.2f%12.2f  %12.2f%12.2f%12.2f   %12.2f%12.2f%3.2f" % tuple(tp) + \
                "%12.2f%12.2f%12.2f" % (chi2_1,chi2_2,chi2_3)
    if fprint is not None:
        print >>fprint,"%12.2f%12.2f%3.2f" % tuple(pars) + \
                "%12.2f%12.2f%12.2f  %12.2f%12.2f%3.2f" % tuple(dp) + \
                "%12.2f%12.2f%12.2f  %12.2f%12.2f%12.2f   %12.2f%12.2f%3.2f" % tuple(tp) +\
                "%12.2f%12.2f%12.2f" % (chi2_1,chi2_2,chi2_3)

    title("Spectrum at %s %s" % (ratos(self.xtora(i)),dectos(self.ytodec(j))) ) 
    legend(loc='best')

def mapplot(plane,cube,vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x):
    
    figure(0)
    clf()
    imshow(plane)

    sp = SpecPlotter(cube,vconv=vconv,xtora=xtora,ytodec=ytodec)
    connect('button_press_event',sp)
 
def splat(filename,vmin=None,vmax=None,button=1,dobaseline=False,exclude=None,smooth=None,order=1,savepre=None,**kwargs):
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
        specplot[:,0,0] = convolve(specplot[:,0,0],hanning(smooth)/hanning(smooth).sum(),'same')

    sp = SpecPlotter(specplot,vconv=vconv,xtora=xtora,ytodec=ytodec)

    sp.plotspec(0,0,button=button,ivconv=ivconv,dv=dv,**kwargs)

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


