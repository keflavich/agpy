# gaussfitter.py
# created by Adam Ginsburg (adam.ginsburg@colorado.edu or keflavich@gmail.com) 3/17/08)
import numpy
from numpy.ma import median
from scipy import optimize,stats,pi
from mpfit import mpfit

"""
To do:
    -turn into a class instead of a collection of objects
    -implement WCS-based gaussian fitting with correct coordinates
"""

def moments(data,circle,rotate,vheight,estimator=median,**kwargs):
    """Returns (height, amplitude, x, y, width_x, width_y, rotation angle)
    the gaussian parameters of a 2D distribution by calculating its
    moments.  Depending on the input parameters, will only output 
    a subset of the above.
    
    If using masked arrays, pass estimator=numpy.ma.median
    """
    total = numpy.abs(data).sum()
    Y, X = numpy.indices(data.shape) # python convention: reverse x,y numpy.indices
    x = (X*numpy.abs(data)).sum()/total
    y = (Y*numpy.abs(data)).sum()/total
    col = data[int(y),:]
    # FIRST moment, not second!
    width_x = numpy.sqrt(numpy.abs((numpy.arange(col.size)-y)*col).sum()/numpy.abs(col).sum())
    row = data[:, int(x)]
    width_y = numpy.sqrt(numpy.abs((numpy.arange(row.size)-x)*row).sum()/numpy.abs(row).sum())
    width = ( width_x + width_y ) / 2.
    height = estimator(data.ravel())
    amplitude = data.max()-height
    mylist = [amplitude,x,y]
    if numpy.isnan(width_y) or numpy.isnan(width_x) or numpy.isnan(height) or numpy.isnan(amplitude):
        raise ValueError("something is nan")
    if vheight==1:
        mylist = [height] + mylist
    if circle==0:
        mylist = mylist + [width_x,width_y]
        if rotate==1:
            mylist = mylist + [0.] #rotation "moment" is just zero...
            # also, circles don't rotate.
    else:  
        mylist = mylist + [width]
    return mylist

def twodgaussian(inpars, circle=0, rotate=1, vheight=1, shape=None):
    """Returns a 2d gaussian function of the form:
        x' = numpy.cos(rota) * x - numpy.sin(rota) * y
        y' = numpy.sin(rota) * x + numpy.cos(rota) * y
        (rota should be in degrees)
        g = b + a * numpy.exp ( - ( ((x-center_x)/width_x)**2 +
        ((y-center_y)/width_y)**2 ) / 2 )

        inpars = [b,a,center_x,center_y,width_x,width_y,rota]
                 (b is background height, a is peak amplitude)

        where x and y are the input parameters of the returned function,
        and all other parameters are specified by this function

        However, the above values are passed by list.  The list should be:
        inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)

        You can choose to ignore / neglect some of the above input parameters 
            unumpy.sing the following options:
            circle=0 - default is an elliptical gaussian (different x, y
                widths), but can reduce the input by one parameter if it's a
                circular gaussian
            rotate=1 - default allows rotation of the gaussian ellipse.  Can
                remove last parameter by setting rotate=0
            vheight=1 - default allows a variable height-above-zero, i.e. an
                additive constant for the Gaussian function.  Can remove first
                parameter by setting this to 0
            shape=None - if shape is set (to a 2-parameter list) then returns
                an image with the gaussian defined by inpars
        """
    inpars_old = inpars
    inpars = list(inpars)
    if vheight == 1:
        height = inpars.pop(0)
        height = float(height)
    else:
        height = float(0)
    amplitude, center_y, center_x = inpars.pop(0),inpars.pop(0),inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    if circle == 1:
        width = inpars.pop(0)
        width_x = float(width)
        width_y = float(width)
        rotate = 0
    else:
        width_x, width_y = inpars.pop(0),inpars.pop(0)
        width_x = float(width_x)
        width_y = float(width_y)
    if rotate == 1:
        rota = inpars.pop(0)
        rota = pi/180. * float(rota)
        rcen_x = center_x * numpy.cos(rota) - center_y * numpy.sin(rota)
        rcen_y = center_x * numpy.sin(rota) + center_y * numpy.cos(rota)
    else:
        rcen_x = center_x
        rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + \
                " circle=%d, rotate=%d, vheight=%d" % (circle,rotate,vheight) )
            
    def rotgauss(x,y):
        if rotate==1:
            xp = x * numpy.cos(rota) - y * numpy.sin(rota)
            yp = x * numpy.sin(rota) + y * numpy.cos(rota)
        else:
            xp = x
            yp = y
        g = height+amplitude*numpy.exp(
            -(((rcen_x-xp)/width_x)**2+
            ((rcen_y-yp)/width_y)**2)/2.)
        return g
    if shape is not None:
        return rotgauss(*numpy.indices(shape))
    else:
        return rotgauss

def gaussfit(data,err=None,params=[],autoderiv=1,return_all=0,circle=0,
        fixed=numpy.repeat(False,7),limitedmin=[False,False,False,False,True,True,True],
        limitedmax=[False,False,False,False,False,False,True],
        usemoment=numpy.array([],dtype='bool'),
        minpars=numpy.repeat(0,7),maxpars=[0,0,0,0,0,0,360],
        rotate=1,vheight=1,quiet=True,returnmp=False,**kwargs):
    """
    Gaussian fitter with the ability to fit a variety of different forms of
    2-dimensional gaussian.
    
    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)
            if not input, these will be determined from the moments of the system, 
            assuming no rotation
        autoderiv=1 - use the autoderiv provided in the lmder.f function (the
            alternative is to us an analytic derivative with lmdif.f: this method
            is less robust)
        return_all=0 - Default is to return only the Gaussian parameters.  See
            below for detail on output
        circle=0 - default is an elliptical gaussian (different x, y widths),
            but can reduce the input by one parameter if it's a circular gaussian
        rotate=1 - default allows rotation of the gaussian ellipse.  Can remove
            last parameter by setting rotate=0.  numpy.expects angle in DEGREES
        vheight=1 - default allows a variable height-above-zero, i.e. an
            additive constant for the Gaussian function.  Can remove first
            parameter by setting this to 0
        usemoment - can choose which parameters to use a moment estimation for.
            Other parameters will be taken from params.  Needs to be a boolean
            array.

    Output:
        Default output is a set of Gaussian parameters with the same shape as
            the input parameters

        Can also output the covariance matrix, 'infodict' that contains a lot
            more detail about the fit (see scipy.optimize.leastsq), and a message
            from leastsq telling what the exit status of the fitting routine was

        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    usemoment=numpy.array(usemoment,dtype='bool')
    params=numpy.array(params,dtype='float')
    if usemoment.any() and len(params)==len(usemoment):
        moment = numpy.array(moments(data,circle,rotate,vheight,**kwargs),dtype='float')
        params[usemoment] = moment[usemoment]
    elif params == [] or len(params)==0:
        params = (moments(data,circle,rotate,vheight,**kwargs))
    if vheight==0:
        vheight=1
        params = numpy.concatenate([[0],params])
        fixed[0] = 1


    # mpfit will fail if it is given a start parameter outside the allowed range:
    for i in xrange(len(params)): 
        if params[i] > maxpars[i] and limitedmax[i]: params[i] = maxpars[i]
        if params[i] < minpars[i] and limitedmin[i]: params[i] = minpars[i]

    if err == None:
        errorfunction = lambda p: numpy.ravel((twodgaussian(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data))
    else:
        errorfunction = lambda p: numpy.ravel((twodgaussian(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data)/err)
    def mpfitfun(data,err):
        if err == None:
            def f(p,fjac=None): return [0,numpy.ravel(data-twodgaussian(p,circle,rotate,vheight)\
                    (*numpy.indices(data.shape)))]
        else:
            def f(p,fjac=None): return [0,numpy.ravel((data-twodgaussian(p,circle,rotate,vheight)\
                    (*numpy.indices(data.shape)))/err)]
        return f

                    
    parinfo = [ 
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"XSHIFT",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"YSHIFT",'error':0},
                {'n':4,'value':params[4],'limits':[minpars[4],maxpars[4]],'limited':[limitedmin[4],limitedmax[4]],'fixed':fixed[4],'parname':"XWIDTH",'error':0} ]
    if vheight == 1:
        parinfo.insert(0,{'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0})
    if circle == 0:
        parinfo.append({'n':5,'value':params[5],'limits':[minpars[5],maxpars[5]],'limited':[limitedmin[5],limitedmax[5]],'fixed':fixed[5],'parname':"YWIDTH",'error':0})
        if rotate == 1:
            parinfo.append({'n':6,'value':params[6],'limits':[minpars[6],maxpars[6]],'limited':[limitedmin[6],limitedmax[6]],'fixed':fixed[6],'parname':"ROTATION",'error':0})

    if autoderiv == 0:
        # the analytic derivative, while not terribly difficult, is less
        # efficient and useful.  I only bothered putting it here because I was
        # instructed to do so for a class project - please ask if you would
        # like this feature implemented
        raise ValueError("I'm sorry, I haven't implemented this feature yet.")
    else:
#        p, cov, infodict, errmsg, success = optimize.leastsq(errorfunction,\
#                params, full_output=1)
        mp = mpfit(mpfitfun(data,err),parinfo=parinfo,quiet=quiet)
    if returnmp:
        return mp
    elif return_all == 0:
        return mp.params
    elif return_all == 1:
        return mp.params,mp.perror

def onedgaussian(x,H,A,dx,w):
    """
    Returns a 1-dimensional gaussian of form
    H+A*numpy.exp(-(x-dx)**2/(2*w**2))
    """
    return H+A*numpy.exp(-(x-dx)**2/(2*w**2))

def onedgaussfit(xax,data,err=None,params=[0,1,0,1],fixed=[False,False,False,False],limitedmin=[False,False,False,True],
        limitedmax=[False,False,False,False],minpars=[0,0,0,0],maxpars=[0,0,0,0],
        quiet=True,shh=True):
    """
    Inputs:
       xax - x axis
       data - y axis
       err - error corresponding to data

       params - Fit parameters: Height of background, Amplitude, Shift, Width
       fixed - Is parameter fixed?
       limitedmin/minpars - set lower limits on each parameter (default: width>0)
       limitedmax/maxpars - set upper limits on each parameter
       quiet - should MPFIT output each iteration?
       shh - output final parameters?

    Returns:
       Fit parameters
       Model
       Fit errors
       chi2
    """

    def mpfitfun(x,y,err):
        if err == None:
            def f(p,fjac=None): return [0,(y-onedgaussian(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-onedgaussian(x,*p))/err]
        return f

    if xax == None:
        xax = numpy.arange(len(data))

    parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0} ,
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"SHIFT",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"WIDTH",'error':0}]

    mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror
    chi2 = mp.fnorm

    if not shh:
        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            print parinfo[i]['parname'],p," +/- ",mpperr[i]
        print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

    return mpp,onedgaussian(xax,*mpp),mpperr,chi2



