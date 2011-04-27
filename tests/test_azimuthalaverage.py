from agpy import azimuthalAverage
from pylab import *


yy,xx = indices([100,100])

rr1 = hypot(xx-50,yy-50)
rr2 = hypot(xx-49.5,yy-49.5)
rr3 = hypot(xx-49.43,yy-49.53)

exp1 = exp(-(rr1**2)/(2.0*50**2))
exp2 = exp(-(rr2**2)/(2.0*50**2))
exp3 = exp(-(rr3**2)/(2.0*50**2))
exp1 /= exp1.max()
exp2 /= exp2.max()
exp3 /= exp3.max()

azr1,azav1 = azimuthalAverage(exp1,center=[50,50],binsize=1.0,returnradii=True)
azr2,azav2 = azimuthalAverage(exp2,center=[49.5,49.5],binsize=1.0,returnradii=True)
azr3,azav3 = azimuthalAverage(exp3,center=[49.43,49.53],binsize=1.0,returnradii=True)
azr1b,azav1b = azimuthalAverage(exp1,center=[50,50],binsize=0.5,returnradii=True)
azr2b,azav2b = azimuthalAverage(exp2,center=[49.5,49.5],binsize=0.5,returnradii=True)
azr3b,azav3b = azimuthalAverage(exp3,center=[49.43,49.53],binsize=0.5,returnradii=True)

subplot(231)
plot(azr1,azav1)
title("Center 50,50, binsize 1")
subplot(232)
plot(azr1b,azav1b)
title("Center 50,50, binsize 0.5")
subplot(233)
plot(azr2,azav2)
title("Center 49.5,49.5, binsize 1")
subplot(234)
plot(azr2b,azav2b)
title("Center 49.5,49.5, binsize 0.5")
subplot(235)
plot(azr3,azav3)
title("Center 49.43,49.53, binsize 1")
subplot(236)
plot(azr3b,azav3b)
title("Center 49.43,49.53, binsize 0.5")

savefig("azimuthalaverage_test.png")
