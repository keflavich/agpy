import numpy as np
import pylab
import matplotlib
old_errsettings = np.geterr()
import pymc # pymc breaks np error settings
np.seterr(**old_errsettings)

def find_percentile(data, pctile):
    sorted_data = np.sort(data.ravel())
    accum_data = sorted_data.cumsum()
    pctiles = accum_data / accum_data.max() * 100.
    return sorted_data[np.argmin(np.abs(pctiles-pctile))]

def errellipse(MC, varname1, varname2, ax=None):
    N = pymc.NormApprox(MC)
    N.fit()
    E =  matplotlib.patches.Ellipse(N.mu[N.__dict__[varname1],
        N.__dict__[varname2]], N.C[N.__dict__[varname1]],
        N.C[N.__dict__[varname2]], 
        (N.C[N.__dict__[varname1], N.__dict__[varname2]][0,1] /
            N.C[N.__dict__[varname1]] * 90.)[0],
        facecolor='none', edgecolor='black')
    if ax is None:
        ax=pylab.gca()
    ax.add_artist(E)

def hist2d(MC, varname1, varname2, varslice=None,
        percentiles=[0.0027,0.0455,0.3173,0.5,0.75],
        colors=[(0.4,0.4,1,0.2),(1,0.4,1,0.5),(1,0.2,0.2,0.5),(0.7,0.1,0.1,1),(0.5,0.05,0.05,1),(0.4,0.05,0.05,0.5)],
        ticklabels=['3$\\sigma$','2$\\sigma$','1$\\sigma$','50%','25%'],
        fignum=1,
        contourcmd=pylab.contourf,
        clear=False,
        colorbar=True,
        doerrellipse=True,
        **kwargs):
    """
    Create a 2D histogram of the MCMC data over some Trace range
    """
    try:
        histvals,xvals,yvals = pylab.histogram2d(MC[varname1].squeeze(),MC[varname2].squeeze(),**kwargs)
    except TypeError:
        if varslice is None:
            histvals,xvals,yvals = pylab.histogram2d(MC.trace(varname1)[:].squeeze(),
                                                     MC.trace(varname2)[:].squeeze(),
                                                     **kwargs)
        else:
            histvals,xvals,yvals = pylab.histogram2d(MC.trace(varname1)[slice(*varslice)].squeeze(),
                                                     MC.trace(varname2)[slice(*varslice)].squeeze(),
                                                     **kwargs)

    levels = [find_percentile(histvals, p*100) for p in percentiles]
    
    pylab.figure(fignum); 
    if clear: pylab.clf(); 

    xax = np.linspace(xvals.min(),xvals.max(),histvals.shape[1])
    yax = np.linspace(yvals.min(),yvals.max(),histvals.shape[0])
    contourcmd(xax, yax, histvals.swapaxes(0,1), levels+[histvals.max()], colors=colors)
    pylab.xlabel(varname1); 
    pylab.ylabel(varname2); 
    if colorbar: 
        cb = pylab.colorbar(); 
        cb.ax.set_yticks(levels); 
        cb.ax.set_yticklabels(ticklabels)

    if doerrellipse:
        errellipse(MC,varname1,varname2)


def gkde_contours(MC, varname1, varname2, varslice=None,
        percentiles=[0.0027,0.0455,0.3173,0.5,0.75],
        colors=[(0.4,0.4,1,0.2),(1,0.4,1,0.5),(1,0.2,0.2,0.75),(1,0.1,0.1,1),(0.8,0.0,0.0,1),(0,0,0,1)],
        ticklabels=['3$\\sigma$','2$\\sigma$','1$\\sigma$','50%','75%'],
        fignum=1,
        ngridpts=101,
        clear=False,):
    """
    """
    import scipy.stats
    data1 = MC.trace(varname1)[slice(*varslice)]
    data2 = MC.trace(varname2)[slice(*varslice)] 
    gkde = scipy.stats.gaussian_kde([data1,data2])
    xvals = np.linspace(data1.min(),data1.max(),ngridpts)
    yvals = np.linspace(data2.min(),data2.max(),ngridpts)
    xx,yy = np.meshgrid(xvals, yvals)

    zz = np.array(gkde.evaluate([xx.flatten(),yy.flatten()])).reshape(xx.shape)

    pylab.figure(fignum); 
    if clear: pylab.clf(); 

    contour(xx, yy, zz, linewidths=1, alpha=.5, cmap=cm.Greys)
    
    pylab.xlabel(varname1); 
    pylab.ylabel(varname2); 
