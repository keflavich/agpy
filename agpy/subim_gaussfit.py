#!/opt/local/bin/ipython
#"""
#Usage:
#./subim_gaussfit.py @file_list.txt x_cen y_cen [outfile.txt]
#%run subim_gaussfit.py @file_list.txt x_cen y_cen [outfile.txt] in an interactive session
# or %run subim_gaussfit.py @file_list.txt @coord_list.txt [outfile.txt]
#"""
import sys
#sys.path.append('/Users/adam/classes/probstat')  #the gaussfitter.py file is in this directory
try:
    import pyfits
except ImportError:
    print "subim_gaussfit requires pyfits"
from numpy import *
from scipy import *
from pylab import *
for k,v in pylab.__dict__.iteritems():  
    if hasattr(v,'__module__'):
        if v.__module__ is None:
            locals()[k].__module__ = 'pylab'
from gaussfitter import *

# read the input files
if len(sys.argv) > 2:
    if sys.argv[1][0] == "@":
        filename = sys.argv[1].strip("@")
        filelist = open(filename).readlines()
    else:
        filelist = [sys.argv[1]]
    if sys.argv[2][0] == "@":
        coord_filename = sys.argv[2].strip("@")
        coord_file = open(coord_filename)
        x_cen = []
        y_cen = []
        for myline in coord_file.readlines():
            x_cen.append(myline.split()[0])
            y_cen.append(myline.split()[1])
        
    elif sys.argv[2][0] != "@":
        x_cen,y_cen = sys.argv[2:4]
    if sys.argv[2][0] == "@" and len(sys.argv) > 3:
        outfile = sys.argv[3]
    elif sys.argv[2][0] != "@" and len(sys.argv) > 4: 
        outfile = sys.argv[4]
else:
    raise ValueError("Wrong number of input parameters.  Input should be of form: " +\
            "\n./subim_gaussfit.py file_list.txt x_cen y_cen ")

# none of these are necessary, they should probably be removed
dates=[]
width=[]
amp=[]
back=[]
ImArr=[]
ModArr=[]
flux=[]

# fit multiple stars function (assumes filelist, x_cen, y_cen are all lists)
# return is a list of lists of lists:
# outermost list the list of stars
#   each star has a list of each time point
#       each time point has a list of parameters
def fitstars(filelist,x_cen,y_cen):
    outdata = []
    if type(x_cen) == type([]):
        for i in xrange(len(x_cen)):
            xc,yc = float(x_cen[i]),float(y_cen[i])
            star_data=[]
            if type(filelist)==type([]):
                for filename in filelist:
                    filename = filename.rstrip('\n').split()[0]
                    if len(filename.split()) > 1:
                        errfilename = filename.split()[1]
                        star_data.append(fitstar(filename,x_cen,y_cen,err=errfilename))
                    else:
                        star_data.append(fitstar(filename,xc,yc))
            else:
                star_data.append(fitstar(filename,xc,yc))
            outdata.append(star_data)
    else:
        star_data=[]
        x_cen,y_cen = float(x_cen),float(y_cen)
        if type(filelist)==type([]):
            for filename in filelist:
                filename = filename.rstrip('\n')
                if len(filename.split()) > 1:
                    errfilename = filename.split()[1]
                    star_data.append(fitstar(filename,x_cen,y_cen,err=errfilename))
                else:
                    star_data.append(fitstar(filename,x_cen,y_cen))
        else:
            star_data.append(fitstar(filename,x_cen,y_cen))
        outdata.append(star_data)
    return outdata


# prints gaussian fit parameters, julian date, and measured flux to file
def printfits(outfilename,filelist,xcen,ycen):
    file = open(outfilename,'w')
    data = fitstars(filelist,xcen,ycen)
    pickle.dump(data,open('.'.join(outfile.split(".")[:-1] + ['pickle.txt']),'w'))
    for i in xrange(len(data)):
        star = data[i]
        print >>file, "# Star %d" % i
        print >>file, "# %15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s " % ('JD','flux','err','height','amplitude','xcen','ycen','xwidth','ywidth','rotation','modelresid','reducedchi2')
        for myline in star:
            print >>file, "%17.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f%15.5f " % tuple(myline)
        print >>file,"\n"

def fitstar(filename,x_cen,y_cen,errname=[]):
    subim_size = 10
    dx = subim_size/2

    file = pyfits.open(filename)
    data = file[0].data
    head = file[0].header
#    date = head['DATE-OBS']
    JD = head['MJD-OBS']
#    angle = int(head['OBJANGLE'])
#    subim = data[x_cen-dx:x_cen+dx,y_cen-dx:y_cen+dx] #star 1
    subim = data[y_cen-dx:y_cen+dx,x_cen-dx:x_cen+dx] #star 1
    b,a,x0,y0,wx,wy,rota = moments(subim,0,1,1)

    if errname == []:
        noise = sqrt(abs(subim))
        radial_weight =  ( ((indices(subim.shape)[0]-x0)**2) + ( (indices(subim.shape)[1]-y0)**2 ) ) ** .5
        errim = noise*radial_weight
    else:
        errfile = pyfits.open(errname)
        errdata = errfile[0].data
        errim = errdata[y_cen-dx:y_cen+dx,x_cen-dx:x_cen+dx]

    parms , cov , infodict , errmsg = gaussfit(subim,err=errim,return_all=1)
    if parms[1] > 2.*subim.max() or parms[1] < subim.max()/2. or parms[0] < subim.mean()/4.:
        ptemp = gaussfit(subim)
        parms , cov , infodict , errmsg = gaussfit(subim,err=errim,params=ptemp,return_all=1)
    b,a,x0,y0,wx,wy,rota = parms
    model = twodgaussian(parms,0,1,1)(*indices(subim.shape))
    chi2 = (((subim - model)/errim)**2).sum()
#    print "date: %s angle: %d height: %.1f amplitude: %.1f  x0: %.1f  y0: %.1f  wx: %.1f wy: %.1f rota: %.1f  chi2: %.1f" % (date,angle,b,a,x0,y0,wx,wy,rota,sum_residuals2)
    wx,wy = abs(wx),abs(wy)
    maskarr = ( asarray (((indices(subim.shape)[0]-x0)**2)<(wx*2) , dtype='int' ) * asarray ( ( (indices(subim.shape)[1]-y0)**2 ) < (wy*2) ,dtype='int'))
    flux = ((maskarr*subim).sum())
    my_error = ((maskarr*errim).sum())
    model_resid = abs(subim-model).sum()
    returnval = tuple([JD]+[flux]+[my_error]+parms.tolist()+[model_resid]+[chi2])
#    return JD,flux,parms,sum_residuals2
    return returnval


printfits(outfile,filelist,x_cen,y_cen)

ParmArr = fitstars(filelist,x_cen,y_cen)
figure(1); clf()
title('FWHM vs Julian Date')
plot(dates,width,'x')
figure(2); clf();
plot(dates,amp,'o')
plot(dates,back,'d')
figure(4); clf();
title('Flux vs. Amplitude')
plot(flux,amp,'x')

