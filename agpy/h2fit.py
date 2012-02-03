from pylab import *
import pylab
for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'
from numpy import *
import agpy
from agpy import readcol
try:
    import pyfits
except ImportError:
    print "pyfits required for h2fit.py"
import timeit
import re
# use pylab's copy instead import copy
from mpfit import mpfit
import gaussfitter

# define physical constants to high precision
h=6.626068e-27
c=2.99792e10
k=1.3806503e-16
e=4.803e-12

tablepath = agpy.__path__[0]+"/h2fit/"

def h2level_energy(V,J):
    """ Returns the theoretical level energy as a function of the
    vibrational (V) and rotational (J) state of the molecule. 
    in units of ergs
    
    Constants are from NIST: 
    http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1000#Diatomic
    (see the bottom of the table)
    """

    We=4401.21
    Be=60.853
    WeXe=121.33 
    De=.0471
    Ae=3.062
    re=.74144

    return h * c * (We*(V+0.5) + Be*(J*(J+1)) - WeXe*(V+.5)**2 - De*J**2*(J+1)**2 - Ae*(V+.5)*(J+1)*J)

# read in rest energies before calling function
resten = readcol(tablepath+'dalgarno1984_table5.txt',verbose=0)

def restwl(vu,vl,ju,jl,calc=False):
    """ Uses energy levels measured by Dabrowski & Herzberg, Can J. Physics, 62,1639,1984 
    vu,vl - upper and lower vibrational states
    ju,jl - upper and lower rotational states 
    returns wavelength in microns
    online versions of this table:
    http://www.astronomy.ohio-state.edu/~depoy/research/observing/molhyd.htm
    http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
    """
    if calc:
        return 1e4*h*c / (h2level_energy(vu,ju) - h2level_energy(vl,jl))
    else:
        if ju >= resten.shape[0] or vu >= resten.shape[1]:
            return 0
        dl = .01/(resten[ju][vu]-resten[jl][vl])
        return dl * 1e6

def linename_to_restwl(linelistfile = tablepath+'linelist.txt',outfile=tablepath+'newlinelist.txt'):

    lines = readcol(linelistfile,fsep='|',twod=False,comment="#")
    outf = open(outfile,'w')

    for line in transpose(lines):
        name = line[0]
        jre = re.compile('\(([0-9]*)\)').search(name)
        if jre == None:
            print >>outf, "%10s|%10s" % (line[0],line[1])
            continue
        else:
            jl = int( jre.groups()[0] )
        if name[4] == 'S':
            ju = jl + 2
        elif name[4] == 'Q': 
            ju = jl
        elif name[4] == 'O':
            ju = jl - 2
        else:
            print >>outf, "%10s|%10s" % (line[0],line[1])
            continue
        vu = int( name[0] )
        vl = int( name[2] )
        rwl = restwl(vu,vl,ju,jl)
        if rwl == 0:
            rwl = float(line[1])
        print >>outf,"%10s|%10.8f" % (name,rwl)



def aval(v,ju,jl):
    """
    Lookup table for Einstein-A value as a function of 
    vibrational level, upper/lower J level
    Values from: http://www.jach.hawaii.edu/UKIRT/astronomy/calib/spec_cal/h2_s.html
    """
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

# atran = pyfits.open('/Users/adam/observations/triplespec/Spextool2/data/atran2000.fits')
atran = readcol(tablepath+'atran.txt')
atran_wl = atran[:,0]*1e4
atran_tr = atran[:,1]
atran_arc = readcol(tablepath+'atran_arcturus.txt')
ARCSORT = argsort(atran_arc[:,0])
atran_arcwl = atran_arc[ARCSORT,0]*1e4
atran_arctr = atran_arc[ARCSORT,1]
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

def readspexspec(image,
    linelistfile = '/Users/adam/work/IRAS05358/code/linelist.txt',
    path_obs='', #'/Users/adam/observations/IRAS05358/UT090108/',
    nameregex='2-1 S\(1\)|1-0 S\([1379028]\)|1-0 Q\([1234]\)|3-2 S\([35]\)|4-3 S\(5\)',
    vlsrcorr=0,
    backsub=False,
    **kwargs):

    regex = re.compile(nameregex)
    im = pyfits.open(path_obs+image)
    wlA = im[0].data[0,:] * 1e4
    data = im[0].data[1,:]

    countsperflux = 2.25e18
    errspec = im[0].data[2,:]

    lines = readcol(linelistfile,fsep='|',twod=False)
    lines[1] = asarray(lines[1],dtype='float')*1e4 # convert microns to angstroms

    specsegments=[]
    for line in transpose(lines):
        wl = float(line[1])
        name = line[0]
        if wl > wlA.min() and wl < wlA.max() and regex.search(name) != None:
            closest = argmin(abs(wlA-wl))
            minind = closest-7
            maxind = closest+7

            savedata = data[minind:maxind]
            if backsub:
                savedata -= median(savedata[savedata<median(savedata)])
            specsegments.append({
                'name':name,
                'linewl':wl,
                'index':closest,
                'minind':minind,
                'maxind':maxind,
                'vlsrcorr':vlsrcorr,
                'wavelength':wlA[minind:maxind],
                'data':savedata,
                'noback':data[minind:maxind]*0,
                'err':errspec[minind:maxind],
                'model':data[minind:maxind]*0
                })

    print "Done finding lines"            

    return specsegments

def readspec(image,noiseimage,
    linelistfile = '/Users/adam/work/IRAS05358/code/linelist.txt',
    path_obs='', #'/Users/adam/observations/IRAS05358/UT090108/',
    noiseaperture=[0,10],
    aperture=[],
    nameregex='2-1 S\(1\)|1-0 S\([1379028]\)|1-0 Q\([1234]\)|3-2 S\([35]\)|4-3 S\(5\)',
    apname='',
    vlsrcorr=0):

    regex = re.compile(nameregex)
    im = pyfits.open(path_obs+image)
    noiseim = pyfits.open(path_obs+noiseimage)
    wlA = im[0].header['CRVAL1'] + im[0].header['CD1_1'] * ( arange(im[0].data.shape[1]) - im[0].header['CRPIX1'] + 1)
    data = im[0].data
    atmospec = median(noiseim[0].data[noiseaperture[0]:noiseaperture[1]],axis=0)

    countsperflux = 2.25e18
    stdatmo = noiseim[0].data[noiseaperture[0]:noiseaperture[1]].std(axis=0)            # std. dev. of non-backsubtracted data
    poisserr = sqrt(abs(noiseim[0].data).mean(axis=0) * countsperflux) / countsperflux  # Poisson noise (very approximate correction)
    errspec = sqrt( stdatmo**2 + 2*poisserr**2 )                 # poisson statistics - once for estimation of the noise, once for the subtraction
    errspec /= atmotrans_vect(wlA)**2                                                  # Weight by inverse of atmospheric transmission^2

    lines = readcol(linelistfile,fsep='|',twod=False)
    lines[1] = asarray(lines[1],dtype='float')*1e4 # convert microns to angstroms

    specsegments=[]
    for line in transpose(lines):
        wl = float(line[1])
        name = line[0]
        if wl > wlA.min() and wl < wlA.max() and regex.search(name) != None:
            closest = argmin(abs(wlA-wl))
            minind = closest-7
            maxind = closest+7

            if len(aperture) != 0:
                savedata = data[aperture[0]:aperture[1],minind:maxind].sum(axis=0)
            else:
                savedata = data[:,minind:maxind]
            specsegments.append({
                'name':name,
                'apname':apname,
                'linewl':wl,
                'index':closest,
                'minind':minind,
                'maxind':maxind,
                'vlsrcorr':vlsrcorr,
                'wavelength':wlA[minind:maxind],
                'data':savedata,
                'noback':atmospec[minind:maxind],
                'err':errspec[minind:maxind],
                'model':data[0,minind:maxind]*0
#                'smoothdata':convolve( data[:,minind:maxind] , hanning(3) , 'same')
                })

    print "Done finding lines"            

    return specsegments

def readstarspec(image,noiseimage,
    linelistfile = '/Users/adam/work/IRAS05358/code/linelist.txt',
    path_obs='/Users/adam/work/IRAS05358/spectra/nearir/', #'/Users/adam/observations/IRAS05358/UT090108/',
    noiseaperture=[0,10],aperture=[]):

    im = pyfits.open(path_obs+image)
    noiseim = pyfits.open(path_obs+noiseimage)
    wlA = im[0].header['CRVAL1'] + im[0].header['CD1_1'] * ( arange(im[0].data.shape[1]) - im[0].header['CRPIX1'] + 1)
    data = im[0].data
    atmospec = median(noiseim[0].data[noiseaperture[0]:noiseaperture[1]],axis=0)
    starspec = data[aperture[0]:aperture[1]].sum(axis=0)

    countsperflux = 2.25e18
    stdatmo = noiseim[0].data[noiseaperture[0]:noiseaperture[1]].std(axis=0)            # std. dev. of non-backsubtracted data
    poisserr = sqrt(abs(noiseim[0].data).mean(axis=0) * countsperflux) / countsperflux  # Poisson noise (very approximate correction)
    errspec = sqrt( stdatmo**2 + 2*poisserr**2 )                 # poisson statistics - once for estimation of the noise, once for the subtraction
    errspec /= atmotrans_vect(wlA)**2                                                  # Weight by inverse of atmospheric transmission^2

    return wlA,starspec,atmospec,errspec

def modelspec(x,T,A,w,dx,op,Ak=0,extinction=False):
    model = x * 0
    A=abs(A)
    # assume width, shift given in velocity:
    w = w*mean(x)/3e5
    dx = dx*mean(x)/3e5
    for v in xrange(1,6):
        for j in xrange(1,14):
            if (j % 2):  # ortho/para multiplier
                mult=op
            else: 
                mult=1
            # S branch
            wl = restwl(v,v-1,j,j-2) * 10**4
            model += A*mult*(2*j+1)*aval(v,j,j-2)*exp(-h2level_energy(v,j)/(k*T)) * exp( - ( x - wl - dx )**2 / (2*w**2) )
            # Q branch
            wl = restwl(v,v-1,j,j) * 10**4
            model += A*mult*(2*(j)+1)*aval(v,j,j)*exp(-h2level_energy(v,j)/(k*T)) * exp( - ( x - wl - dx )**2 / (2*w**2) )
    if extinction:
        # alpha=1.8 comes from Martin & Whittet 1990.  alpha=1.75 from Rieke and Lebofsky 1985
        Al = Ak * (x/22000.0)**(-1.75)
        model *= exp(-Al)
    return model

def testnone(x):
    return int(x != None)

def nonetozero(x):
    if x == None: return 0
    else: return x

def modelpars(fixpars=[None,None,None,None,None],
        minpars=[200.0,0,0,-500.0,1.5],
        maxpars=[15000.0,None,1000.0,500.0,3.5],
        params=[2000,5e-9,31.1,-20.0,3],
        extinction=False,
        **kwargs):

    if len(kwargs) > 0:
        extinction=kwargs['extinction']

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

    if extinction:
        parinfo.append({'n':5,'value':1,'limits':[0,30],'limited':[1,1],'fixed':0,'parname':"AK",'error':0})

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

def modpar(parinfo,fieldnum,value=None,fixed=False,lowlim=None,uplim=None):
    parinfo[fieldnum]['fixed'] = fixed
    if value != None:
        parinfo[fieldnum]['value'] = value
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

def fitspec(specsegments,aperture=None,
        modelspec=modelspec,
        outpars=False,
        vlsrcorr=0,
        extinction=False,
        parinfo=modelpars(),
        quiet=1,
        **kwargs):
    """ fit a model spectrum 
    The model is defined internal to fitspec so that parameters can
    be fixed based on input parameters """

    specsegments = copy(specsegments) # copy.deepcopy(specsegments)
    parinfo = copy(parinfo) # copy.deepcopy(parinfo)

    kwargs = {'extinction':extinction}
    if extinction:
        parinfo=modelpars(extinction=extinction)

    def fitfun(x,y,err):
        return lambda(p): (y-modelspec(x,*p,**kwargs))/err

    def mpfitfun(x,y,err):
        def f(p,fjac=None): return [0,(y-modelspec(x,*p,**kwargs))/err]
        return f

    if aperture is not None:
        miny,maxy = aperture

        fitdata = asarray([ss['data'][miny:maxy,:].sum(axis=0) for ss in specsegments]).ravel()
        fiterr = asarray([ss['err'][:]*sqrt(maxy-miny) for ss in specsegments]).ravel()   # since it is a summed aperture, add in quadrature
    else:
        fitdata = asarray([ss['data'] for ss in specsegments]).ravel()
        fiterr = asarray([ss['err'] for ss in specsegments]).ravel() 

    fitwl = asarray([ss['wavelength'] for ss in specsegments]).ravel()

    try:
        parinfo = parinfo.tolist()
        print "Parinfo was an array.  No idea why, it was never set to one.  Ever."
    except:
        pass

    mp = mpfit(mpfitfun(fitwl,fitdata,fiterr),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror

    for i,p in enumerate(mpp):
        parinfo[i]['value'] = p
        if parinfo[i]['parname'] == 'SHIFT':
            print parinfo[i]['parname'],p+vlsrcorr," +/- ",mpperr[i]
        else:
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

    if aperture is not None:
        for ss in specsegments: ss['data']  = ss['data'][miny:maxy,:].sum(axis=0)
    for ss in specsegments: ss['model'] = modelspec(ss['wavelength'],*mpp)

    if outpars: return specsegments,apfitd,parinfo
    else: return specsegments,apfitd


def plotspec(specsegments,scalefactor=1e17,units="mJy",restwl=-17.5,fignum=0,titletype="line",
        matchaxes=False,vframe='lsr',fitcenter=False):
    nplots = len(specsegments)
    i=1;j=1;fig=figure(fignum);clf();fig.subplots_adjust(left=.05,wspace=.3,hspace=.3,right=.95,top=.95,bottom=.05)
    for ss in specsegments:
        if i % 10 == 0:
            print "Beginning figure %i" % (j)
            i = 1
            draw()
            fig=figure(j); clf();fig.subplots_adjust(left=.05,wspace=.3,hspace=.3,right=.95,top=.95,bottom=.05);
            j+=1
        S=subplot(3,3,i) 
        if titletype=="line":
            title(ss['name']+" "+str(ss['linewl']))
        elif titletype=="ap" or titletype=="aperture":
            title(ss['apname'])
        plotss(ss,scalefactor=scalefactor,units=units,restwl=restwl,vframe=vframe,label=False,fitcenter=fitcenter,atmotrans=False)
        # units: 10^-17 erg/s/cm^2/A
#            axis([xax.min(),xax.max(),segdata.min(),segdata.max()])
        i+=1
    draw()

def plotss(ss,scalefactor=1e17,units="mJy",restwl=-17.5,vframe='lsr',fitcenter=False,
        atmotrans=False,label=True,title=""):

    if atmotrans:
        ax  = axes([.1,.3,.8,.6])
        ax.xaxis.set_visible(False)
        ax2 = axes([.1,.1,.8,.2])
        ax2.plot((atran_wl-ss['linewl'])/ss['linewl']*3e5,atran_tr,'b')
        ax2.plot((atran_arcwl-ss['linewl'])/ss['linewl']*3e5,atran_arctr,'k')
        if label: ax2.set_xlabel("V_LSR (km/s)")
    else: 
        ax=gca()
        if label: ax.set_xlabel("V_LSR (km/s)")
    if len(ss['data'].shape) == 2: segdata = ss['data'].sum(axis=0)
    else: segdata = copy(ss['data'])
    if units=="mJy":
        scalefactor = 1e26 * ss['wavelength']**2  / (c*1e8) * 1.0 # c/l^2 * dl
    elif units=="norm":
        scalefactor = 1./max(segdata)
    segdata *= scalefactor
    xax = (ss['wavelength']-ss['linewl'])/ss['linewl'] * 3e5 # km/s
    if vframe=='lsr': xax += ss['vlsrcorr']
    P1=ax.plot(xax,segdata,linestyle='steps-mid',label='data',color='r')
    if ss['noback'].sum() != 0:
        P2=ax.plot(xax,ss['noback']*scalefactor,linestyle='steps-mid-.',label='raw',color='g')
    PE=ax.errorbar(xax,segdata,yerr=ss['err']*scalefactor,fmt=None,color='r')
    if sum(ss['model']) != 0:
        PM=ax.plot(xax,ss['model']*scalefactor,linestyle='steps-mid--',label='model',color='b')
        PR=ax.plot(xax,segdata-ss['model']*scalefactor,linestyle='steps-mid:',color='k',label='residual')
    xmin,xmax,ymin,ymax=ax.axis()[0:4]; 
    xmin,xmax = xax.min(),xax.max()
    if fitcenter:
        ax.vlines(ss['fitvel']+ss['vlsrcorr'],ymin,ymax,'k',linestyle='-.')
    if restwl != None:
        ax.vlines(restwl,ymin,ymax,'k',linestyle='--')
        ax.axis([xmin,xmax,ymin,ymax])
    if atmotrans: ax2.axis([xmin,xmax,-.05,1.05])
    if label: ax.set_ylabel(units)
    if title: ax.set_title(title)
#    import pdb; pdb.set_trace()

def fitandplotspex(image,
    modelspec=modelspec,
    parinfo=modelpars(),
    nameregex='-',
    vlsrcorr=0,
    extinction=False,
    quiet=1,
    **kwargs):

    if isinstance(image,str):
        ss = readspexspec(image,nameregex=nameregex,**kwargs)
    elif isinstance(image,list):
        ss = image
    ssf,af = fitspec(ss,modelspec=modelspec,parinfo=parinfo,quiet=quiet,vlsrcorr=vlsrcorr,extinction=extinction,**kwargs)
    plotspec(ssf)
    return ssf

def fitandplot(image,noiseimage,
    aperture=[],noiseaperture=[0,10],
    modelspec=modelspec,
    parinfo=modelpars(),
    nameregex='-',
    vlsrcorr=0,
    extinction=False,
    quiet=1):

    ss = readspec(image,noiseimage,noiseaperture=noiseaperture,nameregex=nameregex)
    ssf,af = fitspec(ss,aperture,modelspec=modelspec,parinfo=parinfo,quiet=quiet,vlsrcorr=vlsrcorr,extinction=extinction)
    plotspec(ssf)
    return ssf

def extinctfit(ss): #ss1,ss2):
    ss1 = ss[0]
    ss2 = ss[1]
    trans1 = ss1['name'][4]
    trans2 = ss2['name'][4]
    if trans1 == 'S':
        jupper1 = int(ss1['name'][6])+2
        jupper2 = int(ss2['name'][6])
        vupper1 = int(ss1['name'][0])
        vupper2 = int(ss2['name'][0])
        s = 0
        q = 1
    elif trans1 == 'Q':
        jupper1 = int(ss2['name'][6])+2
        jupper2 = int(ss1['name'][6])
        vupper1 = int(ss1['name'][0])
        vupper2 = int(ss2['name'][0])
        s = 1
        q = 0
    else:
        print "Error: wrong transition type"
        return
    if jupper1 != jupper2 or vupper1!=vupper2:
        print "Error: upper levels not matched"
        return

    ASdivAQ = aval(vupper1,jupper1,jupper1-2) / aval(vupper1,jupper1,jupper1)   # S / Q
    WSdivWQ = restwl(vupper1,vupper1-1,jupper1,jupper1-2) / restwl(vupper1,vupper1-1,jupper1,jupper1)

    midwl = [median(ss1['wavelength']),median(ss2['wavelength'])]
    pguess1 = [0,ss1['data'].max()*1e17,midwl[0],5]
    pguess2 = [0,ss2['data'].max()*1e17,midwl[1],5]
    p1,fit1,perr1,gof1 = gaussfitter.onedgaussfit(ss1['wavelength'],ss1['data']*1e17,ss1['err']*1e17,params=pguess1,minpars=[-1,0,19000,0],limitedmin=[True,True,True,True])
    p2,fit2,perr2,gof2 = gaussfitter.onedgaussfit(ss2['wavelength'],ss2['data']*1e17,ss2['err']*1e17,params=pguess2,minpars=[-1,0,19000,0],limitedmin=[True,True,True,True])

    ss[0]['model'] = fit1*1e-17
    ss[1]['model'] = fit2*1e-17

    if s == 0:   # want S / Q
        SdivQ = fit1.sum()/fit2.sum()  
    else:
        SdivQ = fit2.sum()/fit1.sum()

    dtau12 = -log(SdivQ/ASdivAQ*WSdivWQ)  # see Moore, Lumsden, Ridge, Puxley 2005
                                     # alpha = 1.8 comes from Martin & Whittet 1990
    tau1 = dtau12 / ( (midwl[s]/midwl[q])**-1.8 - 1 )
    AL = 100**(.2) * log10(exp(1)) * tau1 # conversion from optical depth to magnitude
    # AK = AL * lambda^-alpha
    AK = AL * (midwl[q]/22000.0)**-1.8 # again alpha=1.8 comes from Martin & Whittet 1990.  alpha=1.75 from Rieke and Lebofsky 1985

    return AK,AL,SdivQ,ASdivAQ

def linefit(ss):
    pguess = [0,ss['data'].max()*1e17,0,20]
    wl = (ss['wavelength']-ss['linewl'])/ss['linewl']*3e5
    pa,fit,perr,gof = gaussfitter.onedgaussfit(wl,ss['data']*1e17,ss['err']*1e17,params=pguess,
            minpars=[-1,0,-200,0],limitedmin=[True,True,True,True])
#    import pdb; pdb.set_trace()
    ss['model'] = fit*1e-17
    return pa,perr,gof

def fitalllines(specsegments,apname='aperture'):

    fitstruct = []
    print "%10s%20s%20s%20s" % ('name','flux','velocity','width')
    for ss in specsegments:
        pa,perr,gof = linefit(ss)

        goodregion = (ss['model'] - pa[0]*1e-17) > 5e-17

        ss['fitvel'] = pa[2]
        ss['fitvelerr'] = perr[2]
        ss['fitwidth'] = pa[3]
        ss['fitwidtherr'] = perr[3]
        ss['flux'] = ss['model'].sum()
        ss['fluxerr'] = (ss['err'][goodregion]).sum()
        print "%10s%10.4g(%10.4g)%10.4g(%10.4g)%10.4g(%10.4g)" % (
                ss['name'],ss['model'].sum(),(sqrt((ss['err'][goodregion]**2).sum())),pa[2]+ss['vlsrcorr'],perr[2],pa[3],perr[3])
        fitstruct.append({'name':ss['name'],
          'linewl':ss['linewl'],
          'flux':ss['model'].sum(),
          'error':sqrt((ss['err'][goodregion]**2).sum()),
          'velocity':pa[2]+ss['vlsrcorr'],
          'vel_err':perr[2],
          'width':pa[3],
          'width_err':perr[3],
          'goodness':gof,
          'red_chi2':gof/goodregion.sum() } )

    plotspec(specsegments)
    return fitstruct

def extinctandplot(specsegments):
    ssj0 = [ None , None ]
    ssj1 = [ None , None ]
    ssj2 = [ None , None ]
    for ss in specsegments:
        if ss['name'][:8]=='1-0 S(0)':
            ssj0[0] = ss
        elif ss['name'][:8]=='1-0 S(1)':
            ssj1[0] = ss
        elif ss['name'][:8]=='1-0 S(2)':
            ssj2[0] = ss
        elif ss['name'][:8]=='1-0 Q(2)':
            ssj0[1] = ss
        elif ss['name'][:8]=='1-0 Q(3)':
            ssj1[1] = ss
        elif ss['name'][:8]=='1-0 Q(4)':
            ssj2[1] = ss
    if ssj0[0] != None and ssj0[1] != None:
        AK0,AJ0,fr0,ar0 = extinctfit(ssj0)
        plotspec(ssj0,fignum=0)
        print "J=2 fluxratio=%g AK=%g AV=%g" % (fr0,AK0,AK0/.112) # .112 from Rieke and Lebofsky 1985
    if ssj1[0] != None and ssj1[1] != None:
        AK1,AJ1,fr1,ar1 = extinctfit(ssj1)
        plotspec(ssj1,fignum=1)
        print "J=3 fluxratio=%g AK=%g AV=%g" % (fr1,AK1,AK1/.112)
    if ssj2[0] != None and ssj2[1] != None:
        AK2,AJ2,fr2,ar2 = extinctfit(ssj2)
        plotspec(ssj2,fignum=2)
        print "J=4 fluxratio=%g AK=%g AV=%g" % (fr2,AK2,AK2/.112)
    return AK0,AK1,AK2

def printfslist(fslist,apnamelist,outfilename='line_measurements.csv',fsep=","):
    outf= open(outfilename,'w')
    printnames = [",%20s" % line['name'].rstrip() for line in fslist[0]]
    printwls = [",%20g" % (float(line['linewl'])/1e4) for line in fslist[0]]
    print >>outf,"          ",fsep,fsep.join(printnames)
    print >>outf,"          ",fsep,fsep.join(printwls  )
    for i,fs in enumerate(fslist):
        printfluxes = [fsep.join(["%10.2g" % float(line['flux']),"(%8.2g)" % float(line['error'])]) for line in fs]
        print >>outf,"%10s" % apnamelist[i],fsep,fsep.join(printfluxes)
    print >>outf,""
    for i,fs in enumerate(fslist):
        printfluxes = [fsep.join(["%10.2g" % float(line['velocity']),"(%8.2g)" % float(line['vel_err'])]) for line in fs]
        print >>outf,"%10s" % apnamelist[i],fsep,fsep.join(printfluxes)
    print >>outf,""
    for i,fs in enumerate(fslist):
        printfluxes = [fsep.join(["%10.2g" % float(line['goodness']),"%10.2g" % float(line['red_chi2'])]) for line in fs]
        print >>outf,"%10s" % apnamelist[i],fsep,fsep.join(printfluxes)

    outf.close()
    

