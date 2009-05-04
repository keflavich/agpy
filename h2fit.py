from pylab import *
from numpy import *
from readcol import readcol

h=6.626068e-27
c=2.99792e10
k=1.3806503e-16
e=4.803e-12

We=4401.21
Be=60.853
WeXe=121.33 
De=.0471
Ae=3.062
re=.74144

def h2level_energy(V,J):
    return h * c * (We*(V+0.5) + Be*(J*(J+1)) - WeXe*(V+.5)**2 - De*J**2*(J+1)**2 - Ae*(V+.5)*(J+1)*J)

def flux(wl,T,gu,gl):
    return gu / gl * exp(-k*T)

def restwl(vu,ju,jl):
    """ Essentially a lookup table.  
    To do: replace with either a theoretical
    calculation or reading a text file"""
    if vu == 1:
        if ju-jl==2:
            if jl==0 : return  2.2235
            elif jl==1 : return  2.1218
            elif jl==2 : return  2.0338
            elif jl==3 : return  1.9576
            elif jl==4 : return  1.8920
            elif jl==5 : return  1.8358
            elif jl==6 : return  1.7880
            elif jl==7 : return  1.7480
            elif jl==8 : return  1.7147
            elif jl==9 : return  1.6877
            elif jl==10: return  1.6665
            elif jl==11: return  1.6504
            else: return 0
        elif ju-jl==0:
            if jl==1: return  2.4066
            elif jl==2: return  2.4134
            elif jl==3: return  2.4237
            elif jl==4: return  2.4375
            elif jl==5: return  2.4548
            elif jl==6: return  2.4756
            elif jl==7: return  2.5001
            else: return 0
    elif vu == 2:
        if jl==0  and ju-jl==2: return  2.3556
        elif jl==1  and ju-jl==2: return  2.2477
        elif jl==2  and ju-jl==2: return  2.1542
        elif jl==3  and ju-jl==2: return  2.0735
        elif jl==4  and ju-jl==2: return  2.0041
        elif jl==5  and ju-jl==2: return  1.9449
        else: return 0
    elif vu == 3:
        if jl==0  and ju-jl==2: return  2.5014
        elif jl==1  and ju-jl==2: return  2.3864
        elif jl==2  and ju-jl==2: return  2.2870
        elif jl==3  and ju-jl==2: return  2.2014
        elif jl==4  and ju-jl==2: return  2.1280
        elif jl==5  and ju-jl==2: return  2.0656
        elif jl==6  and ju-jl==2: return  2.0130
        elif jl==7  and ju-jl==2: return  1.9692
        else: return 0
    else: return 0

def aval(v,ju,jl):
    if v==1:
        if jl==0 and ju-jl==2: return 2.53e-7 
        elif jl==1 and ju-jl==2: return 3.47e-7 
        elif jl==2 and ju-jl==2: return 3.98e-7 
        elif jl==3 and ju-jl==2: return 4.21e-7 
        elif jl==4 and ju-jl==2: return 4.19e-7 
        elif jl==5 and ju-jl==2: return 3.96e-7 
        elif jl==6 and ju-jl==2: return 3.54e-7 
        elif jl==7 and ju-jl==2: return 2.98e-7 
        elif jl==8 and ju-jl==2: return 2.34e-7 
        elif jl==9 and ju-jl==2: return 1.68e-7 
        elif jl==1 and ju-jl==0: return 4.29e-7 
        elif jl==2 and ju-jl==0: return 3.03e-7 
        elif jl==3 and ju-jl==0: return 2.78e-7 
        elif jl==4 and ju-jl==0: return 2.65e-7 
        else: return 0
    elif v==2:
        if jl==0 and ju-jl==2: return 3.68e-7 
        elif jl==1 and ju-jl==2: return 4.98e-7 
        elif jl==2 and ju-jl==2: return 5.60e-7 
        elif jl==3 and ju-jl==2: return 5.77e-7 
        elif jl==4 and ju-jl==2: return 5.57e-7 
        else: return 0
    elif v==3:
        if jl==0 and ju-jl==2: return 3.88e-7 
        elif jl==1 and ju-jl==2: return 5.14e-7 
        elif jl==2 and ju-jl==2: return 5.63e-7 
        elif jl==3 and ju-jl==2: return 5.63e-7 
        elif jl==4 and ju-jl==2: return 5.22e-7 
        elif jl==5 and ju-jl==2: return 4.50e-7 
        else: return 0
    else: return 0

aval_vect=vectorize(aval)

def relflux(T,g,wl,jupper,jlower,vupper):
    f0 = s_1 # 1.00 # flux of S(0) 1-0
    f0 = 21.0 * exp(-h2level_energy(1,3)/(k*T)) * aval(1,3,1)
    f = g * aval_vect(vupper,jupper,jlower) * exp( -h2level_energy(vupper,jupper) / (k*T)) / f0

    return f

# what am I trying to fit?
# set some 'zero point flux' for, say, 1-0 S(0) 
# flux as a function of [transition] should be something....
# if (Q branch) then zero is Q(1)
# if (S branch) then zero is S(0) depending on which V....


import readcol
import pyfits
import timeit
import re
import copy
import mpfit

atran = pyfits.open('/Users/adam/observations/triplespec/Spextool2/data/atran2000.fits')
atran_wl = atran[0].data[0,:]*1e4
atran_tr = atran[0].data[1,:]
def atmotrans(x):
    """ returns the atmospheric transmission at the given wavelength (in angstroms) """
    closest = argmin(abs(atran_wl-x))
    if atran_wl[closest] < x:
        m = (atran_tr[closest+1]-atran_tr[closest])/(atran_wl[closest+1]-atran_wl[closest])
        b = atran_tr[closest]
        y = m * (x-atran_wl[closest]) + b
    elif atran_wl[closest] > x:
        m = (atran_tr[closest]-atran_tr[closest-1])/(atran_wl[closest]-atran_wl[closest-1])
        b = atran_tr[closest-1]
        y = m * (x-atran_wl[closest-1]) + b
    else:
        y = atran_tr[closest]
    return y

atmotrans_vect = vectorize(atmotrans)

def readspec(image,noiseimage,
    linelistfile = '/Users/adam/work/IRAS05358/code/linelist.txt',
    path_obs='/Users/adam/observations/IRAS05358/UT090108/',
    noiseaperture=[0,10],
    nameregex='-'):

    regex = re.compile(nameregex)
    im = pyfits.open(path_obs+image)
    noiseim = pyfits.open(path_obs+noiseimage)
    wlA = im[0].header['CRVAL1'] + im[0].header['CD1_1'] * ( arange(im[0].data.shape[1]) - im[0].header['CRPIX1'])
    data = im[0].data
    atmospec = median(noiseim[0].data[noiseaperture[0]:noiseaperture[1]],axis=0)

    countsperflux = 2.25e18
    stdatmo = noiseim[0].data[noiseaperture[0]:noiseaperture[1]].std(axis=0)            # std. dev. of non-backsubtracted data
    poisserr = sqrt(abs(noiseim[0].data).mean(axis=0) * countsperflux) / countsperflux  # Poisson noise (very approximate correction)
    errspec = sqrt( stdatmo**2 + 2*poisserr**2 )                 # poisson statistics - once for estimation of the noise, once for the subtraction
    errspec /= atmotrans_vect(wlA)**2                                                  # Weight by inverse of atmospheric transmission^2

    lines = readcol.readcol(linelistfile,fsep='|',twod=False,dtype='S')
    lines[1] = asarray(lines[1],dtype='float')*1e4 # convert microns to angstroms

    specsegments=[]
    for line in transpose(lines):
        wl = float(line[1])
        name = line[0]
        if wl > wlA.min() and wl < wlA.max() and regex.search(name) != None:
            closest = argmin(abs(wlA-wl))
            minind = closest-7
            maxind = closest+7
            specsegments.append({
                'name':name,
                'linewl':wl,
                'index':closest,
                'minind':minind,
                'maxind':maxind,
                'wavelength':wlA[minind:maxind],
                'data':data[:,minind:maxind],
                'noback':atmospec[minind:maxind],
                'err':errspec[minind:maxind],
                'model':data[:,minind:maxind]*0
#                'smoothdata':convolve( data[:,minind:maxind] , hanning(3) , 'same')
                })

    print "Done finding lines"            

    return specsegments


def modelspec(x,T,A,w,dx,op):
    model = x * 0
    A=abs(A)
    # assume width, shift given in velocity:
    w = w*mean(x)/3e5
    dx = dx*mean(x)/3e5
    for v in xrange(1,5):
        for j in xrange(1,16):
            if (j % 2):  # ortho/para multiplier
                mult=op
            else: 
                mult=1
            # S branch
            wl = restwl(v,j,j-2) * 10**4
            model += A*mult*(2*j+1)*aval(v,j,j-2)*exp(-h2level_energy(v,j)/(k*T)) * exp( - ( x - wl - dx )**2 / (2*w**2) )
            # Q branch
            wl = restwl(v,j,j) * 10**4
            model += A*mult*(2*(j)+1)*aval(v,j,j)*exp(-h2level_energy(v,j)/(k*T)) * exp( - ( x - wl - dx )**2 / (2*w**2) )
    return model

def testnone(x):
    return int(x != None)

def nonetozero(x):
    if x == None: return 0
    else: return x

def modelpars(fixpars=[None,None,None,None,None],
        minpars=[200.0,0,0,-500.0,1.5],
        maxpars=[15000.0,None,1000.0,500.0,3.5],
        params=[2000,5e-9,3.1,-2.0,3]):

    limitedmin = map(testnone,minpars)
    limitedmax = map(testnone,maxpars)
    fixed = map(testnone,fixpars)
    minpars = map(nonetozero,minpars)
    maxpars = map(nonetozero,maxpars)
    fixpars = map(nonetozero,fixpars)

    parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"TEMPERATURE",'error':0} ,
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"SCALE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"WIDTH",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"SHIFT",'error':0},
                {'n':4,'value':params[4],'limits':[minpars[4],maxpars[4]],'limited':[limitedmin[4],limitedmax[4]],'fixed':fixed[4],'parname':"ORTHOtoPARA",'error':0}]
    return parinfo

def twoTmodelspec(x,T1,A1,w1,dx1,op1,T2,A2,w2,dx2,op2):
    sumspec = modelspec(x,T1,A1,w1,dx1,op1) + modelspec(x,T2,A2,w2,dx2,op2)
    return sumspec

def twoTmodelpars(fixpars=[None,None,None,None,None,None,None,None,None,None],
        minpars=[200.0,0,0,-500.0,1.0,200.0,0,0,-500.0,1.0],
        maxpars=[15000.0,None,1000.0,500.0,5.0,15000.0,None,1000.0,500.0,5.0],
        params=[2000,5e-9,3.1,-2.0,3,2000,5e-9,3.1,-2.0,3]):

    limitedmin = map(testnone,minpars)
    limitedmax = map(testnone,maxpars)
    fixed = map(testnone,fixpars)
    minpars = map(nonetozero,minpars)
    maxpars = map(nonetozero,maxpars)
    fixpars = map(nonetozero,fixpars)

    parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"TEMPERATURE",'error':0} ,
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"SCALE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"WIDTH",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"SHIFT",'error':0},
                {'n':4,'value':params[4],'limits':[minpars[4],maxpars[4]],'limited':[limitedmin[4],limitedmax[4]],'fixed':fixed[4],'parname':"ORTHOtoPARA",'error':0},
                {'n':5,'value':params[5],'limits':[minpars[5],maxpars[5]],'limited':[limitedmin[5],limitedmax[5]],'fixed':fixed[5],'parname':"TEMPERATURE",'error':0},
                {'n':6,'value':params[6],'limits':[minpars[6],maxpars[6]],'limited':[limitedmin[6],limitedmax[6]],'fixed':fixed[6],'parname':"SCALE",'error':0},
                {'n':7,'value':params[7],'limits':[minpars[7],maxpars[7]],'limited':[limitedmin[7],limitedmax[7]],'fixed':fixed[7],'parname':"WIDTH",'error':0},
                {'n':8,'value':params[8],'limits':[minpars[8],maxpars[8]],'limited':[limitedmin[8],limitedmax[8]],'fixed':fixed[8],'parname':"SHIFT",'error':0},
                {'n':9,'value':params[9],'limits':[minpars[9],maxpars[9]],'limited':[limitedmin[9],limitedmax[9]],'fixed':fixed[9],'parname':"ORTHOtoPARA",'error':0}]
    return parinfo

def modpar(parinfo,fieldnum,value,fixed=False,lowlim=None,uplim=None):
    parinfo[fieldnum]['value'] = value
    parinfo[fieldnum]['fixed'] = fixed
    if lowlim != None:
        parinfo[fieldnum]['limits'][0] = lowlim
        parinfo[fieldnum]['limited'][0] = True
    if uplim != None:
        parinfo[fieldnum]['limits'][1] = uplim
        parinfo[fieldnum]['limited'][1] = True

def showpars(parinfo):
    print "%20s  %12s%12s%10s%10s%10s" % ("Parameter","Value","Error","Low","High","Frozen?")
    for i,p in enumerate(parinfo):
        print "%4i%16s  %12.6g%12.6g%10.4g%10.4g%10i" % (i,p['parname'],p['value'],p['error'],p['limits'][0],p['limits'][1],p['fixed'])

def fitspec(specsegments,aperture,
        modelspec=modelspec,
        parinfo=modelpars(),
        outpars=False,
        quiet=1):
    """ fit a model spectrum 
    The model is defined internal to fitspec so that parameters can
    be fixed based on input parameters """

    specsegments = copy.deepcopy(specsegments)
    parinfo = copy.deepcopy(parinfo)

    def fitfun(x,y,err):
        return lambda(p): (y-modelspec(x,*p))/err

    def mpfitfun(x,y,err):
        def f(p,fjac=None): return [0,(y-modelspec(x,*p))/err]
        return f

    miny,maxy = aperture

    fitwl = asarray([ss['wavelength'] for ss in specsegments]).ravel()
    fitdata = asarray([ss['data'][miny:maxy,:].sum(axis=0) for ss in specsegments]).ravel()
    fiterr = asarray([ss['err'][:]*sqrt(maxy-miny) for ss in specsegments]).ravel()   # since it is a summed aperture, add in quadrature

    mp = mpfit.mpfit(mpfitfun(fitwl,fitdata,fiterr),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror

    for i,p in enumerate(mpp):
        parinfo[i]['value'] = p
        print parinfo[i]['parname'],p," +/- ",mpperr[i]
#    print "Temperature: ",mpp[0]," Shift: ",mpp[3],mpp[3]*mean(fitwl)/3e5," Width: ",mpp[2],mpp[2]*mean(fitwl)/3e5,"   Ortho/Para: ",mpp[4]
#    print "ERRORS: Temperature: ",mpperr[0]," Shift: ",mpperr[3],mpperr[3]*mean(fitwl)/3e5," Width: ",mpperr[2],mpperr[2]*mean(fitwl)/3e5,"   Ortho/Para: ",mpperr[4]
    print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(fitdata)," DOF:",len(fitdata)-len(mpp)

    apfitd = {
        'params':mpp,
        'parerr':mpperr,
        'parinfo':parinfo,
        'wl':fitwl,
        'data':fitdata,
        'model':modelspec(fitwl,*mpp),
        'err':fiterr
        }

    for ss in specsegments: ss['data']  = ss['data'][miny:maxy,:].sum(axis=0)
    for ss in specsegments: ss['model'] = modelspec(ss['wavelength'],*mpp)

    if outpars: return specsegments,apfitd,parinfo
    else: return specsegments,apfitd


def plotspec(specsegments,scalefactor=1e17):
    nplots = len(specsegments)
    i=1;j=1;figure(0);clf();
    for ss in specsegments:
        if i % 10 == 0:
            print "Beginning figure %i" % (j)
            i = 1
            draw()
            figure(j); clf();
            j+=1
        S=subplot(3,3,i) 
        title(ss['name']+" "+str(ss['linewl']))
        # units: 10^-17 erg/s/cm^2/A
        segdata = ss['data']*scalefactor
        xax = (ss['wavelength']-ss['linewl'])/ss['linewl'] * 3e5 # km/s
        P1=plot(xax,segdata,linestyle='steps-mid',label='data',color='r')
        P2=plot(xax,ss['noback']*scalefactor,linestyle='steps-mid-.',label='raw',color='g')
        PE=errorbar(xax,segdata,yerr=ss['err']*scalefactor,fmt=None,color='r')
        PM=plot(xax,ss['model']*scalefactor,linestyle='steps-mid--',label='model',color='b')
        PR=plot(xax,segdata-ss['model']*scalefactor,linestyle='steps-mid:',color='k',label='residual')
#            axis([xax.min(),xax.max(),segdata.min(),segdata.max()])
        i+=1
    draw()

def fitandplot(image,noiseimage,
    apertures=[],noiseaperture=[0,10],
    modelspec=modelspec,
    parinfo=modelpars(),
    nameregex='-'):

    ss = readspec(image,noiseimage,noiseaperture=noiseaperture,nameregex=nameregex)
    ssf,af = fitspec(ss,apertures,modelspec=modelspec,parinfo=parinfo)
    plotspec(ssf)


show()        
