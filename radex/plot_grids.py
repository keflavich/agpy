from pylab import *
import pyfits
from agpy import readcol,asinh_norm
import matplotlib
import sys

def plot_radex(filename,ngridpts=100,ncontours=50,plottype='ratio',transition="noname",thirdvarname="Temperature",
    cutnumber=None,cutvalue=10,vmin=None,vmax=None,logscale=False,**kwargs):

    names,props = readcol(filename,twod=False,names=True)
    temperature,density,column,tex1,tex2,tau1,tau2,ratio = props


    if thirdvarname == "Temperature":
      firstvar = density
      secondvar = column
      thirdvar = temperature
      if cutnumber is not None:
        cutvalue = unique(thirdvar)[int(cutnumber)]
      firstlabel = "log$(n_{H_2}) ($cm$^{-3})$"
      secondlabel = "log$(N_{H_2CO}) ($cm$^{-2})$"
      savetype = "DenCol_T=%iK" % cutvalue
      graphtitle = "T = %g K" % cutvalue
    elif thirdvarname == "Density":
      firstvar = temperature
      secondvar = column
      thirdvar = density
      if cutnumber is not None:
        cutvalue = unique(thirdvar)[int(cutnumber)]
      firstlabel = "Temperature (K)"
      secondlabel = "log$(N_{H_2CO}) ($cm$^{-2})$"
      savetype = "TemCol_n=1e%gpercc" % cutvalue
      graphtitle = "n = %g cm$^{-3}$" % (10**cutvalue)
    elif thirdvarname == "Column":
      secondvar = density
      firstvar = temperature
      thirdvar = column
      if cutnumber is not None:
        cutvalue = unique(thirdvar)[int(cutnumber)]
      secondlabel = "log$(n_{H_2}) ($cm$^{-3})$"
      firstlabel = "Temperature (K)"
      savetype = "TemDen_N=1e%gpersc" % cutvalue
      graphtitle = "N = %g cm$^{-2}$" % (10**cutvalue)

    if plottype == 'ratio':
      cblabel = "$F_{1-1} / F_{2-2}$"
    elif plottype == 'tau1':
      cblabel = "$\\tau_{1-1}$"
    elif plottype == 'tau2':
      cblabel = "$\\tau_{2-2}$"
    elif plottype == 'tex1':
      cblabel = "$\\T_{ex}(1-1)$"
    elif plottype == 'tex2':
      cblabel = "$\\T_{ex}(2-2)$"

    varfilter = (thirdvar==cutvalue)
    if varfilter.sum() == 0:
      raise ValueError("Cut value %g does not match any of %s values" % (cutvalue, thirdvarname))

    nx = len(unique(firstvar))
    ny = len(unique(secondvar))
    if firstvar is temperature:
      firstarr = logspace(log10(firstvar.min()),log10(firstvar.max()),nx)
    else:
      firstarr = linspace(firstvar.min(),firstvar.max(),nx)
    secondarr = linspace(secondvar.min(),secondvar.max(),ny)

    exec('plotdata = %s' % plottype)

    plot_grid = griddata(firstvar[varfilter],secondvar[varfilter],plotdata[varfilter],firstarr,secondarr)
    
    if vmax:
      plot_grid[plot_grid > vmax] = vmax
    if vmin:
      plot_grid[plot_grid > vmin] = vmin
    if logscale:
      plot_grid = log10(plot_grid)

    figure(1)
    clf()
    conlevs = logspace(-3,1,ncontours)
    contourf(firstarr,secondarr,plot_grid,conlevs,norm=matplotlib.colors.LogNorm())#,**kwargs) #,norm=asinh_norm.AsinhNorm(**kwargs),**kwargs)
    xlabel(firstlabel)
    ylabel(secondlabel)
    title(graphtitle)
    cb = colorbar()
    cb.set_label(cblabel)
    cb.set_ticks([1e-3,1e-2,1e-1,1,1e1])
    cb.set_ticklabels([1e-3,1e-2,1e-1,1,1e1])
    savefig("%s_%s_%s.png" % (savetype,plottype,transition))

def gridcube(filename,outfilename,plottype='tau1',transition='noname'):

    names,props = readcol(filename,twod=False,names=True)
    temperature,density,column,tex1,tex2,tau1,tau2,ratio = props

    exec('plotdata = %s' % plottype)

    nx = len(unique(density))
    ny = len(unique(column))
    nz = len(unique(temperature))

    densarr = linspace(density.min(),density.max(),nx)
    colarr  = linspace(column.min(),column.max(),ny)

    newarr = zeros([nz,ny,nx])

    for iT,T in enumerate(unique(temperature)):
      varfilter = temperature==T
      newarr[iT,:,:] = griddata(density[varfilter],column[varfilter],plotdata[varfilter],densarr,colarr)

    newfile = pyfits.PrimaryHDU(newarr)
    newfile.header.update('CRVAL3' ,  log10(min(temperature)) )
    newfile.header.update('CRPIX3' ,  1 )
    newfile.header.update('CTYPE3' ,  'LOG-TEMP' )
    newfile.header.update('CD3_3' , log10(unique(temperature)[1]) - log10(unique(temperature)[0]) )
    newfile.header.update('CRVAL1' ,  min(densarr) )
    newfile.header.update('CRPIX1' ,  1 )
    newfile.header.update('CD1_1' , densarr[1]-densarr[0] )
    newfile.header.update('CTYPE1' ,  'LOG-DENS' )
    newfile.header.update('CRVAL2' ,  min(colarr) )
    newfile.header.update('CRPIX2' ,  1 )
    newfile.header.update('CD2_2' , colarr[1]-colarr[0] )
    newfile.header.update('CTYPE2' ,  'LOG-COLU' )
    newfile.writeto(outfilename,clobber=True)


if __name__ == "__main__": 

    filename = sys.argv[1]
    if len(sys.argv) > 2:
        transition = sys.argv[2]
    else:
        transition = filename[:7]

    if len(sys.argv) > 3:
        plottype = sys.argv[3]
    else:
        plottype = 'ratio'

    if len(sys.argv) > 4:
        cutnumber = sys.argv[4]
    else:
        cutnumber = 0

    plot_radex(filename,transition=transition,plottype=plottype,cutnumber=cutnumber,thirdvarname="Temperature")
    plot_radex(filename,transition=transition,plottype=plottype,cutnumber=cutnumber,thirdvarname="Density")
    plot_radex(filename,transition=transition,plottype=plottype,cutnumber=cutnumber,thirdvarname="Column")

    show()
