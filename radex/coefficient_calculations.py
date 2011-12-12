#!/usr/bin/env python
# data files come from ftp://cdsarc.u-strasbg.fr/pub/cats/J/A+A/493/687/
from agpy import readcol
import numpy

h2co_oo = readcol('h2co_oo.dat',asStruct=True,comment='%',skipline=3)
h2co_op = readcol('h2co_op.dat',asStruct=True,comment='%',skipline=3)
levels = readcol('o-h2co_levels.dat',asStruct=True)

nlevelso = len(h2co_oo.Jl)
h2co_oo.__dict__['llevelnum'] = numpy.zeros(nlevelso)
h2co_oo.__dict__['ulevelnum'] = numpy.zeros(nlevelso)
nlevelsp = len(h2co_op.Jl)
h2co_op.__dict__['llevelnum'] = numpy.zeros(nlevelsp)
h2co_op.__dict__['ulevelnum'] = numpy.zeros(nlevelsp)

temperatures = numpy.arange(5,105,5)

def R(a0,a1,a2,a3,a4,T):
    """
    Troscompt et al (2009) coefficients using Faure et al (2004) equation:
    log10(R) = sum(a_n T^{-n/6})
    where n=0..4, R is presumably cm^3 s^-1
    """
    return a0 + a1*T**(-1./6.) + a2*T**(-2./6.) + a3*T**(-3./6.) + a4*T**(-4./6.)

outf = open('o-h2co_troscompt.dat','w')

print >>outf, "!MOLECULE"
print >>outf, "o-H2CO"
print >>outf, "!MOLECULAR WEIGHT"
print >>outf, "30.0"
print >>outf, "!NUMBER OF ENERGY LEVELS"
print >>outf, "10"
print >>outf, "!LEVEL + ENERGIES(cm^-1) + WEIGHT + J_Kp_Ko"
print >>outf, "    1       10.539000      3.0      1_1_1"
print >>outf, "    2       10.700100      3.0      1_1_0"
print >>outf, "    3       15.236900      5.0      2_1_2"
print >>outf, "    4       15.720200      5.0      2_1_1"
print >>outf, "    5       22.282200      7.0      3_1_3"
print >>outf, "    6       23.248700      7.0      3_1_2"
print >>outf, "    7       31.672900      9.0      4_1_4"
print >>outf, "    8       33.283500      9.0      4_1_3"
print >>outf, "    9       43.406700     11.0      5_1_5"
print >>outf, "   10       45.822000     11.0      5_1_4"
print >>outf, "!NUMBER OF RADIATIVE TRANSITIONS"
print >>outf, "13"
print >>outf, "!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)"
print >>outf, "  1    2    1      3.564e-09          4.829660     15.4"
print >>outf, "  2    3    1      5.304e-05        140.839502     21.9"
print >>outf, "  3    4    3      3.208e-08         14.488479     22.6"
print >>outf, "  4    4    2      6.472e-05        150.498334     22.6"
print >>outf, "  5    5    3      2.271e-04        211.211468     32.1"
print >>outf, "  6    6    5      1.283e-07         28.974805     33.4"
print >>outf, "  7    6    4      2.772e-04        225.697775     33.4"
print >>outf, "  8    7    5      5.883e-04        281.526929     45.6"
print >>outf, "  9    8    7      3.564e-07         48.284547     47.9"
print >>outf, " 10    8    6      7.178e-04        300.836635     47.9"
print >>outf, " 11    9    7      1.202e-03        351.768645     62.5"
print >>outf, " 12   10    9      8.017e-07         72.409090     65.9"
print >>outf, " 13   10    8      1.466e-03        375.893216     65.9"
print >>outf, "!NUMBER OF COLL PARTNERS"
print >>outf, "2"
# % Lines (13+NLEV+NLIN) - (14+NLEV+NLIN): collision partner ID and reference. Valid identifications are: 1=H2, 2=para-H2, 3=ortho-H2, 4=electrons, 5=H, 6=He, 7=H+. 
print >>outf, "!COLLISIONS BETWEEN"
print >>outf, "3 o-H2CO - o-H2 from Troscompt (2009)"
print >>outf, "!NUMBER OF COLL TRANS"
print >>outf, "45"
print >>outf, "!NUMBER OF COLL TEMPS"
print >>outf, "%i" % (len(temperatures))
print >>outf, "!COLL TEMPS"
print >>outf, " "+" ".join(["%6.1f" % T for T in temperatures])
print >>outf, "!TRANS+ UP+ LOW+ COLLRATES(cm^3 s^-1)"

for ii in xrange(nlevelso):
    llevelname = "%i_%i_%i" % (h2co_oo.Jl[ii],h2co_oo.Kal[ii],h2co_oo.Kcl[ii])
    ulevelname = "%i_%i_%i" % (h2co_oo.Ju[ii],h2co_oo.Kau[ii],h2co_oo.Kcu[ii])
    h2co_oo.llevelnum[ii] = levels.LEVEL[levels.J_Kp_Ko == llevelname]
    h2co_oo.ulevelnum[ii] = levels.LEVEL[levels.J_Kp_Ko == ulevelname]
    collrates = [R(h2co_oo.a0[ii],h2co_oo.a1[ii],h2co_oo.a2[ii],h2co_oo.a3[ii],h2co_oo.a4[ii],T) for T in temperatures]
    print >>outf, "%4i %4i %4i" % (ii+1,h2co_oo.ulevelnum[ii],h2co_oo.llevelnum[ii]),
    print >>outf, " ".join(["%10.3e" % (10**C) for C in collrates])

# % Lines (13+NLEV+NLIN) - (14+NLEV+NLIN): collision partner ID and reference. Valid identifications are: 1=H2, 2=para-H2, 3=ortho-H2, 4=electrons, 5=H, 6=He, 7=H+. 
print >>outf, "!COLLISIONS BETWEEN"
print >>outf, "2 o-H2CO - p-H2 from Troscompt (2009)"
print >>outf, "!NUMBER OF COLL TRANS"
print >>outf, "45"
print >>outf, "!NUMBER OF COLL TEMPS"
print >>outf, "%i" % (len(temperatures))
print >>outf, "!COLL TEMPS"
print >>outf, " ".join(["%0.1f" % T for T in temperatures])
print >>outf, "!TRANS+ UP+ LOW+ COLLRATES(cm^3 s^-1)"

for ii in xrange(nlevelsp):
    llevelname = "%i_%i_%i" % (h2co_op.Jl[ii],h2co_op.Kal[ii],h2co_op.Kcl[ii])
    ulevelname = "%i_%i_%i" % (h2co_op.Ju[ii],h2co_op.Kau[ii],h2co_op.Kcu[ii])
    h2co_op.llevelnum[ii] = levels.LEVEL[levels.J_Kp_Ko == llevelname]
    h2co_op.ulevelnum[ii] = levels.LEVEL[levels.J_Kp_Ko == ulevelname]
    collrates = [R(h2co_op.a0[ii],h2co_op.a1[ii],h2co_op.a2[ii],h2co_op.a3[ii],h2co_op.a4[ii],T) for T in temperatures]
    print >>outf, "%i %i %i " % (ii+1,h2co_op.ulevelnum[ii],h2co_op.llevelnum[ii]),
    print >>outf, " ".join(["%9.3e" % (10**C) for C in collrates])

outf.close()
