import numpy as np

def powerlaw_sample(prob, alpha, minval, maxval):
    """
    """
    # PDF = X_0 x^-alpha
    # CDF = integral(X_0 x^-alpha,xmin,x) = [X_0 / (-alpha+1) x^(-alpha+1)]_xmin^xmax
    #   1 = X_0/(1-alpha)*(xmax^(1-alpha)-xmin^(1-alpha))
    X0 =  (1.-alpha)/(maxval**(1.-alpha)-minval**(1.-alpha))
    #cdf = X0/(1.-alpha) * (x**(1.-alpha)-minval**(1.-alpha))
    #return cdf
    x = ((1.-prob)*(1.-alpha)/X0 + minval**(1.-alpha))**(1./(1.-alpha))
    return x

def ellipses(N, walpha=2, halpha=2, wmin=1., wmax=512., hmin=1., hmax=512.,
        maxx=512., maxy=512.):
    """
    Create N ellipses from powerlaw sampling with the specified parameters
    """

    widths  = powerlaw_sample(np.random.rand(N), walpha, wmin, wmax)
    heights = powerlaw_sample(np.random.rand(N), halpha, hmin, hmax)
    position_angles = np.random.rand(N)*180.
    xcenters = np.random.rand(N)*maxx
    ycenters = np.random.rand(N)*maxy
    
    return zip(xcenters,ycenters,widths,heights,position_angles)


if __name__=="__main__":
    
    data = np.zeros([512,512])
    yy,xx = np.indices(data.shape)

    N=3500.
    for (x,y,w,h,p),amp in zip(ellipses(N,walpha=3),powerlaw_sample(np.random.rand(N),2,0.01,1)):
        th = p/180.*np.pi
        a = np.cos(th)**2/(2.*w**2) + np.sin(th)**2/(2.*h**2)
        b = -np.sin(2*th)/(4.*w**2) + np.sin(2*th)/(4.*h**2)
        c = np.sin(th)**2/(2.*w**2) + np.cos(th)**2/(2.*h**2)
        #print x,y,w,h,p
        data += amp * np.exp(-(xx-x)**2*a-b*(xx-x)*(yy-y)-(yy-y)**2*c)

    import pylab
    pylab.clf()
    pylab.imshow(data)



