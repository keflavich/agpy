import numpy as np
import pylab
import matplotlib
old_errsettings = np.geterr()
try:
    import pymc # pymc breaks np error settings
except ImportError:
    pass
np.seterr(**old_errsettings)

def find_percentile(data, pctile):
    sorted_data = np.sort(data.ravel())
    accum_data = sorted_data.cumsum()
    pctiles = accum_data / accum_data.max() * 100.
    return sorted_data[np.argmin(np.abs(pctiles-pctile))]

def errellipse(MC, varname1, varname2, ax=None):
    N = pymc.NormApprox(MC)
    N.fit()
    E = matplotlib.patches.Ellipse(N.mu[N.__dict__[varname1],
                                   N.__dict__[varname2]],
                                   N.C[N.__dict__[varname1]],
                                   N.C[N.__dict__[varname2]],
                                   (N.C[N.__dict__[varname1],
                                    N.__dict__[varname2]][0,1] /
                                    N.C[N.__dict__[varname1]] * 90.)[0],
                                   facecolor='none',
                                   edgecolor='black')
    if ax is None:
        ax=pylab.gca()
    ax.add_artist(E)

def hist2d(MC, varname1, varname2, varslice=None,
           percentiles=[0.0027,0.0455,0.3173,0.5,0.75],
           colors=[(0.4,0.4,1,0.2),(1,0.4,1,0.5),(1,0.2,0.2,0.5),(0.7,0.1,0.1,1),(0.5,0.05,0.05,1),(0.4,0.05,0.05,0.5)],
           ticklabels=['3$\\sigma$','2$\\sigma$','1$\\sigma$','50%','25%'],
           axis=None,
           fignum=1,
           contourcmd=pylab.contourf,
           clear=False,
           colorbar=True,
           doerrellipse=False,
           chain=None,
           **kwargs):
    """
    Create a 2D histogram of the MCMC data over some Trace range
    """
    try: # if input is just a dict of arrays
        if varslice is None:
            histvals,xvals,yvals = pylab.histogram2d(MC[varname1].squeeze(),MC[varname2].squeeze(),**kwargs)
        else:
            histvals,xvals,yvals = pylab.histogram2d(MC[varname1][slice(*varslice)].squeeze(),MC[varname2][slice(*varslice)].squeeze(),**kwargs)
    except TypeError:
        if varslice is None:
            histvals,xvals,yvals = pylab.histogram2d(
                MC.trace(varname1,chain=chain)[:].squeeze(),
                MC.trace(varname2,chain=chain)[:].squeeze(),
                **kwargs)
        else:
            histvals,xvals,yvals = pylab.histogram2d(
                MC.trace(varname1,chain=chain)[slice(*varslice)].squeeze(),
                MC.trace(varname2,chain=chain)[slice(*varslice)].squeeze(),
                **kwargs)

    levels = [find_percentile(histvals, p*100) for p in percentiles]
    
    if axis is None:
        pylab.figure(fignum)
        if clear:
            pylab.clf()
        axis = pylab.gca()

    xax = np.linspace(xvals.min(),xvals.max(),histvals.shape[1])
    yax = np.linspace(yvals.min(),yvals.max(),histvals.shape[0])
    if axis is not None:
        contourcmd = eval('axis.'+contourcmd.__name__)
    cntr = contourcmd(xax, yax, histvals.swapaxes(0,1), levels+[histvals.max()], colors=colors)
    # hack to fix opacity
    axis.set_xlabel(varname1)
    axis.set_ylabel(varname2)
    if colorbar:
        try:
            cb = pylab.colorbar(cntr, ax=axis)
            cb.ax.set_yticks(levels)
            cb.ax.set_yticklabels(ticklabels)
        except Exception as e:
            print "Colorbar failed with exception ",e

    if doerrellipse:
        errellipse(MC,varname1,varname2)
    return axis


def gkde_contours(MC, varname1, varname2, varslice=None,
                  percentiles=[0.0027,0.0455,0.3173,0.5,0.75],
                  colors=[(0.4,0.4,1,0.2),(1,0.4,1,0.5),(1,0.2,0.2,0.75),(1,0.1,0.1,1),(0.8,0.0,0.0,1),(0,0,0,1)],
                  ticklabels=['3$\\sigma$','2$\\sigma$','1$\\sigma$','50%','75%'],
                  fignum=1,
                  ngridpts=101,
                  clear=False,):
    """
    Contours for kernel densit estimate... to compare to real contours
    """
    import scipy.stats
    data1 = MC.trace(varname1)[slice(*varslice)]
    data2 = MC.trace(varname2)[slice(*varslice)]
    gkde = scipy.stats.gaussian_kde([data1,data2])
    xvals = np.linspace(data1.min(),data1.max(),ngridpts)
    yvals = np.linspace(data2.min(),data2.max(),ngridpts)
    xx,yy = np.meshgrid(xvals, yvals)

    zz = np.array(gkde.evaluate([xx.flatten(),yy.flatten()])).reshape(xx.shape)

    pylab.figure(fignum)
    if clear:
        pylab.clf()

    pylab.contour(xx, yy, zz, linewidths=1, alpha=.5, cmap=matplotlib.cm.Greys)
    
    pylab.xlabel(varname1)
    pylab.ylabel(varname2)


def plot_mc_hist(MC, field, varslice=None, onesided=True, bins=50, chain=None,
        axis=None, lolim=False, legloc='best', legend=True, **kwargs):
    """
    Plot a histogram with 1,2,3-sigma bars
    """
    try:
        field_data = MC[field].squeeze()
    except TypeError:
        field_data = MC.trace(field,chain=chain)[:]
    if varslice is not None:
        field_data = field_data[slice(*varslice)]

    field_stats = {'mean': field_data.mean()}
    if onesided:
        #field_stats = MC.trace(field,chain=chain).stats(quantiles=[68.2689,95.44997,99.7300,50])
        quantiles = {1:68.2689,2:95.44997,3:99.7300,'m':50}
        if lolim:
            quantiles = {k:100-q for k,q in quantiles.iteritems()}
        field_stats['quantiles'] = {k:np.percentile(field_data,q) for k,q in quantiles.iteritems()}
    else:
        #field_stats = MC.trace(field,chain=chain).stats(quantiles=[0.135,2.275,15.866,84.134,97.725,99.865,50])
        field_stats['quantiles'] = {q:np.percentile(field_data,q) for q in [0.135,2.275,15.866,84.134,97.725,99.865,50]}

    vpts = field_stats['quantiles']
    if axis is None:
        ax = pylab.gca()
    else:
        ax = axis
    #field_data_sorted = np.sort(field_data)
    h,l,p = ax.hist(field_data,bins=bins,histtype='stepfilled',**kwargs)
    if kwargs.get('normed'):
        ylim = [0,h.max()*1.01]
    else:
        ylim = ax.get_ylim()
    #fieldlen = len(field_data)
    if onesided:
        ax.vlines(vpts[1], *ylim,linewidth=3, alpha=0.5, color='k',label="$1\\sigma$")
        ax.vlines(vpts[2],*ylim,linewidth=3, alpha=0.5, color='r',label="$2\\sigma$")
        ax.vlines(vpts[3], *ylim,linewidth=3, alpha=0.5, color='g',label="$3\\sigma$")
    else:
        ax.vlines(field_stats['mean'],*ylim,color='k', linestyle='--', linewidth=3, alpha=0.5, label="$\mu$")
        ax.vlines(vpts[50],*ylim, color='b', linestyle='--', linewidth=3, alpha=0.5, label="$\mu_{1/2}$")
        ax.vlines([vpts[15.866],vpts[84.134]],*ylim,color='k',linewidth=3, alpha=0.5, label="$1\\sigma$")
        ax.vlines([vpts[02.275],vpts[97.725]],*ylim,color='r',linewidth=3, alpha=0.5, label="$2\\sigma$")
        ax.vlines([vpts[00.135],vpts[99.865]],*ylim,color='g',linewidth=3, alpha=0.5, label="$3\\sigma$")
    ax.set_ylim(*ylim)
    if legend:
        ax.legend(loc=legloc)

    return ax

def autocorr_diagnostics(mc):
    traces = mc.db._traces
    ntraces = len(traces)
    npanels = np.ceil(np.sqrt(ntraces))

    for ii,(k,v) in enumerate(traces.iteritems()):
        if v[:].ndim > 1: 
            d = v[:,0].squeeze()
        else:
            d = v[:].squeeze()
        pylab.subplot(npanels, npanels, ii+1)
        ft = np.fft.fft(d)
        ac = np.fft.ifft(ft*ft[::-1])
        frq = np.fft.fftfreq(ac.size)
        pylab.plot(frq,ac,',')

def trace_diagnostics(mc):
    traces = mc.db._traces
    ntraces = len(traces)
    npanels = np.ceil(np.sqrt(ntraces))

    for ii,(k,v) in enumerate(traces.iteritems()):
        if v[:].ndim > 1: 
            d = v[:,0].squeeze()
        else:
            d = v[:].squeeze()
        pylab.subplot(npanels, npanels, ii+1)
        pylab.plot(d,',')
