import numpy

def downsample(myarr,factor):
    xs,ys = myarr.shape
    crarr = myarr[:xs-(xs % int(factor)),:ys-(ys % int(factor))]
    dsarr = numpy.concatenate([[crarr[i::factor,j::factor] 
        for i in range(factor)] 
        for j in range(factor)]).mean(axis=0)
    return dsarr
