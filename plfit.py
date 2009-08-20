# intended to implement a power-law fitting routine as specified in.....
# that paper Erik sent me in July 2009.  
# http://www.santafe.edu/~aaronc/powerlaws/
#
# The MLE for the power-law alpha is very easy to derive given knowledge
# of the lowest value at which a power law holds, but that point is 
# difficult to derive and must be acquired iteratively.

"""
plfit.py
Requires pylab (matplotlib), which requires numpy

example use:
from plfit import plfit

MyPL = plfit(mydata)
MyPL.plotpdf(log=True)


"""

from pylab import *
import numpy

class plfit:
    """
    A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
    from http://www.santafe.edu/~aaronc/powerlaws/

    See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
    in empirical data" SIAM Review, to appear (2009). (arXiv:0706.1062)
    http://arxiv.org/abs/0706.1062

    The output "alpha" is defined such that p(x) ~ (x/xmin)^-alpha
    """

    def __init__(self,x,**kwargs):
        """
        Initializes and fits the power law.  Can pass "quiet" to turn off 
        output (except for warnings)
        """
        self.data = x
        self.plfit(**kwargs)


    def alpha_(self,x):
        def alpha(xmin,x=x):
            """
            given a sorted data set and a minimum, returns power law MLE fit
            data is passed as a keyword parameter so that it can be vectorized
            """
            x = x[x>=xmin]
            n = sum(x>=xmin)
            a = float(n) / sum(log(x/xmin))
            return a
        return alpha

    def kstest_(self,x):
        def kstest(xmin,x=x):
            """
            given a sorted data set and a minimum, returns power law MLE ks-test w/data
            data is passed as a keyword parameter so that it can be vectorized
            """
            x = x[x>=xmin]
            n = float(len(x))
            a = float(n) / sum(log(x/xmin))
            cx = arange(n,dtype='float')/float(n)
            cf = 1-(xmin/x)**a
            ks = max(abs(cf-cx))
            return ks
        return kstest
    

    # should probably use a decorator here
    def plfit(self,nosmall=False,finite=False,quiet=False,silent=False):
        """
        A Python implementation of the Matlab code http://www.santafe.edu/~aaronc/powerlaws/plfit.m
        from http://www.santafe.edu/~aaronc/powerlaws/

        See A. Clauset, C.R. Shalizi, and M.E.J. Newman, "Power-law distributions
        in empirical data" SIAM Review, to appear (2009). (arXiv:0706.1062)
        http://arxiv.org/abs/0706.1062
        """
        x = self.data
        xmins = unique(x)[1:-1]
        #xmins = xmins[1:-1]
        z  = sort(x)
        av  = asarray( map(self.alpha_(z),xmins) ,dtype='float')
        dat = asarray( map(self.kstest_(z),xmins),dtype='float')
        if nosmall:
            # test to make sure the number of data points is high enough
            # to provide a reasonable s/n on the computed alpha
            sigma = (av-1)/sqrt(arange(len(av),0,-1))
            dat = dat[sigma<0.1]
        D     = min(dat);
        xmin  = xmins[argmin(dat)]
        z     = x[x>=xmin]
        n     = len(z)
        alpha = 1 + n / sum( log(z/xmin) )
        if finite:
            alpha = alpha*(n-1.)/n+1./n
        if n < 50 and not finite and not silent:
            print '(PLFIT) Warning: finite-size bias may be present. n=%i' % n
        ks = max(abs( arange(n)/float(n) - (1-(xmin/z)**alpha) ))
        L = n*log((alpha-1)/xmin) - alpha*sum(log(z/xmin));
        #requires another map... Larr = arange(len(unique(x))) * log((av-1)/unique(x)) - av*sum
        self._likelihood = L
        self._av = av
        self._xmin_kstest = dat
        self._xmin = xmin
        self._alpha= alpha
        self._alphaerr = (alpha-1)/sqrt(n)
        self._ks = ks

        if not quiet:
            print "xmin: %g  n(>xmin): %i  alpha: %g +/- %g  Likelihood: %g  ks: %g" % (xmin,n,alpha,self._alphaerr,L,ks)

        return xmin,alpha


    def plotcdf(self,x=None,xmin=None,alpha=None,**kwargs):
        """
        Plots CDF and powerlaw
        """
        x=self.data
        xmin=self._xmin
        alpha=self._alpha

        x=sort(x)
        n=len(x)
        xcdf = arange(n,0,-1,dtype='float')/float(n)

        q = x[x>=xmin]
        fcdf = (q/xmin)**(1-alpha)
        nc = xcdf[argmax(x>=xmin)]
        fcdf_norm = nc*fcdf

        loglog(x,xcdf,marker='+',color='k',**kwargs)
        loglog(q,fcdf_norm,color='r',**kwargs)

    def plotpdf(self,x=None,xmin=None,alpha=None,nbins=50,dolog=False,dnds=False,**kwargs):
        """
        Plots PDF and powerlaw....
        """
        x=self.data
        xmin=self._xmin
        alpha=self._alpha

        x=sort(x)
        n=len(x)

        gca().set_xscale('log')
        gca().set_yscale('log')

        if dnds:
            hb = histogram(x,bins=logspace(log10(min(x)),log10(max(x)),nbins))
            h = hb[0]
            b = hb[1]
            db = hb[1][1:]-hb[1][:-1]
            h = h/db
            plot(b[:-1],h,drawstyle='steps-post',color='k',**kwargs)
            #alpha -= 1
        elif dolog:
            hb = hist(x,bins=logspace(log10(min(x)),log10(max(x)),nbins),log=True,fill=False,edgecolor='k',**kwargs)
            alpha -= 1
            h,b=hb[0],hb[1]
        else:
            hb = hist(x,bins=linspace((min(x)),(max(x)),nbins),fill=False,edgecolor='k',**kwargs)
            h,b=hb[0],hb[1]
        b = b[1:]

        q = x[x>=xmin]
        px = (alpha-1)/xmin * (q/xmin)**(-alpha)

        arg = argmin(abs(b-xmin))
        plotloc = (b>xmin)*(h>0)
        norm = median( h[plotloc] / ((alpha-1)/xmin * (b[plotloc]/xmin)**(-alpha))  )
        px = px*norm

        loglog(q,px,'r',**kwargs)
        vlines(xmin,0.1,max(px),colors='r',linestyle='dashed')

        gca().set_xlim(min(x),max(x))

    def test_pl(self,niter=1e3,**kwargs):
        """
        Monte-Carlo test to determine whether distribution is consistent with a power law

        Runs through niter iterations of a sample size identical to the input sample size.

        Will randomly select values from the data < xmin.  The number of values selected will
        be chosen from a uniform random distribution with p(<xmin) = n(<xmin)/n.

        Once the sample is created, is fit using above methods, then the best fit is used to
        compute a Kolmogorov-Smirnov statistic.  The KS stat distribution is compared to the 
        KS value for the fit to the actual data, and p = fraction of random ks values greater
        than the data ks value is computed.  If p<.1, the data may be inconsistent with a 
        powerlaw.  A data set of n(>xmin)>100 is required to distinguish a PL from an exponential,
        and n(>xmin)>~300 is required to distinguish a log-normal distribution from a PL.
        For more details, see figure 4.1 and section

        **WARNING** This can take a very long time to run!  Execution time scales as 
        niter * setsize

        """
        xmin = self._xmin
        alpha = self._alpha
        niter = int(niter)

        ntail = sum(self.data >= xmin)
        ntot = len(self.data)
        nnot = ntot-ntail              # n(<xmin)
        pnot = nnot/float(ntot)        # p(<xmin)
        nonpldata = self.data[self.data<xmin]
        nrandnot = sum( rand(ntot) < pnot ) # randomly choose how many to sample from <xmin
        nrandtail = ntot - nrandnot         # and the rest will be sampled from the powerlaw

        ksv = []
        for i in xrange(niter):
            # first, randomly sample from power law
            # with caveat!  
            nonplind = floor(rand(nrandnot)*nnot).astype('int')
            fakenonpl = nonpldata[nonplind]
            randarr = rand(nrandtail)
            fakepl = randarr**(1/(1-alpha)) * xmin 
            fakedata = concatenate([fakenonpl,fakepl])
            # second, fit to powerlaw
            TEST = plfit(fakedata,quiet=True,silent=True,nosmall=True,**kwargs)
            ksv.append(TEST._ks)
        
        ksv = array(ksv)
        p = (ksv>self._ks).sum() / float(niter)
        self._pval = p
        self._ks_rand = ksv

        return p,ksv

def plfit_lsq(x,y):
    """
    Returns A and B in y=Ax^B
    http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
    """
    n = len(x)
    btop = n * (log(x)*log(y)).sum() - (log(x)).sum()*(log(y)).sum()
    bbottom = n*(log(x)**2).sum() - (log(x).sum())**2
    b = btop / bbottom
    a = ( log(y).sum() - b * log(x).sum() ) / n

    A = exp(a)
    return A,b
