# Fit a PSF of type Airy or Gaussian...
from gaussfitter import twodgaussian,moments
import numpy
import scipy
from numpy import pi
from mpfit import mpfit

def _airy_func(rr, amplitude=1.0, width=1.0):
    """
    For a simple radially symmetric airy function, returns the value at a given
    (normalized) radius
    """
    return amplitude * (2.0 * scipy.special.j1(rr/width) / (rr/width))**2

def _gaussian_func(rr, amplitude=1.0, sigma=1.0):
    """
    For a simple radially symmetric Gaussian function, returns the value at a given
    (normalized) radius
    """
    return amplitude * numpy.exp(-(rr**2) / (2.0 * sigma**2) )

def airy(inpars, circle=True, rotate=False, vheight=True, shape=None, fwhm=False):
    """Returns a 2d Airy *function* of the form:
        x' = numpy.cos(rota) * x - numpy.sin(rota) * y
        y' = numpy.sin(rota) * x + numpy.cos(rota) * y
        (rota should be in degrees)
        radius = sqrt( (x'-xcen)^2 + (y'-ycen)^2 )
        g = b + a * 2.0*BesselJ1( radius ) / radius
            (with a correction for the divide-by-zero in the center)

        inpars = [b,a,center_x,center_y,width_x,width_y,rota]
                 (b is background height, a is peak amplitude)

        where x and y are the input parameters of the returned function,
        and all other parameters are specified by this function

        However, the above values are passed by list.  The list should be:
        inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)

        You can choose to ignore / neglect some of the above input parameters 
            unumpy.sing the following options:
            circle=1 - default is a circular Airy Disk.  An elliptical Airy is
                possible, but probably not physically motivated (unless it's
                sampled onto a stretched grid?).
            rotate=0 - default allows rotation of the gaussian ellipse.  Can
                remove last parameter by setting rotate=0
            vheight=1 - default allows a variable height-above-zero, i.e. an
                additive constant for the Airy function.  Can remove first
                parameter by setting this to 0
            shape=None - if shape is set (to a 2-parameter list) then returns
                an image with the gaussian defined by inpars
        fwhm - if set, assumes the Width parameters input are FWHM widths, so
            they'll be converted to "Sigma" widths by s = FWHM/2.0/1.61633
            (http://en.wikipedia.org/wiki/Airy_disk 
            and http://home.fnal.gov/~neilsen/notebook/astroPSF/astroPSF.html)
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

    if fwhm:
        width_x /= 2.0 * 1.61633
        width_y /= 2.0 * 1.61633
            
    def rotairy(x,y):
        if rotate==1:
            xp = x * numpy.cos(rota) - y * numpy.sin(rota)
            yp = x * numpy.sin(rota) + y * numpy.cos(rota)
        else:
            xp = x
            yp = y
        rr = numpy.sqrt(((rcen_x-xp)/width_x)**2+
             ((rcen_y-yp)/width_y)**2)
        # http://en.wikipedia.org/wiki/Airy_disk
        airy_func = (2.0 * scipy.special.j1(rr) / rr)**2
        airy_func[rr==0] = 1.0
        airy = height + amplitude * airy_func

        return airy
    if shape is not None:
        return rotairy(*numpy.indices(shape))
    else:
        return rotairy

def psffit(data,err=None,params=[],autoderiv=True,return_all=False,circle=True,
        fixed=numpy.repeat(False,7),limitedmin=[False,False,False,False,True,True,True],
        limitedmax=[False,False,False,False,False,False,True],
        usemoment=numpy.array([],dtype='bool'),
        minpars=numpy.repeat(0,7),maxpars=[0,0,0,0,0,0,360],
        rotate=0,vheight=1,quiet=True,returnmp=False,
        returnfitimage=False,
        psffunction=airy, 
        extra_pars=None,
        return_parinfo=False,
        **kwargs):
    """
    PSF fitter with the ability to fit a variety of different forms of
    2-dimensional gaussian OR an Airy.
    This code is mostly directly copied from gaussfitter.py and presents 
    yet another argument for me turning this into a class...
    
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
        return_all=0 - Default is to return only the Gaussian parameters.  
                   1 - fit params, fit error
        returnfitimage - returns (best fit params,best fit image)
        returnmp - returns the full mpfit struct
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
        extra_pars - If your psffunction requires extra parameters, pass their
            parinfo dictionaries through this variable

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

    if err is None:
        errorfunction = lambda p: numpy.ravel((psffunction(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data))
    else:
        errorfunction = lambda p: numpy.ravel((psffunction(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data)/err)
    def mpfitfun(data,err):
        if err is None:
            def f(p,fjac=None): return [0,numpy.ravel(data-psffunction(p,circle,rotate,vheight)\
                    (*numpy.indices(data.shape)))]
        else:
            def f(p,fjac=None): return [0,numpy.ravel((data-psffunction(p,circle,rotate,vheight)\
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

    if extra_pars:
        for P in extra_pars:
            parinfo.append(P)

    if autoderiv is False:
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
        returns = (mp)
    elif return_parinfo:
        returns = (parinfo)
    elif return_all is False:
        returns = mp.params
    elif return_all:
        returns = mp.params,mp.perror
    if returnfitimage:
        fitimage = psffunction(mp.params,circle,rotate,vheight)(*numpy.indices(data.shape))
        returns = (returns,fitimage)
    return returns
