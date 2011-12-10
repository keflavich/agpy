from pylab import *
from agpy import correlate2d,psds,convolve
import line_profiler

LP = line_profiler.LineProfiler()

exponent = 9

xx,yy = indices([2**exponent*1.5,2**exponent*1.5])
rr1 = sqrt((xx-2**exponent)**2+(yy-2**exponent)**2)
rr2 = sqrt((xx-2**(exponent-1)*1.5)**2+(yy-2**(exponent-1)*1.5)**2)

ee1 = exp(-rr1**2/2.0**4)
ee2 = exp(-rr2**2/2.0**4)

acorr1pc = correlate2d(ee1,ee1,fft_pad=True)
acorr1nc = correlate2d(ee1,ee1,fft_pad=False)
acorr2pc = correlate2d(ee2,ee2,fft_pad=True)
acorr2nc = correlate2d(ee2,ee2,fft_pad=False)

acorr1pn = correlate2d(ee1,ee1,fft_pad=True,crop=False)
acorr1nn = correlate2d(ee1,ee1,fft_pad=False,crop=False)
acorr2pn = correlate2d(ee2,ee2,fft_pad=True,crop=False)
acorr2nn = correlate2d(ee2,ee2,fft_pad=False,crop=False)

LP.add_function(convolve)
aconv1n = LP.runcall(convolve,ee1,ee1,fft_pad=False)
aconv2n = LP.runcall(convolve,ee2,ee2,fft_pad=False)
print "STATS FOR fft_pad=False"
LP.print_stats()
aconv1p = LP.runcall(convolve,ee1,ee1,fft_pad=True)
aconv2p = LP.runcall(convolve,ee2,ee2,fft_pad=True)
print "STATS FOR fft_pad=True"
LP.print_stats()

xcorr1pc = correlate2d(ee1,ee2,fft_pad=True)
xcorr1nc = correlate2d(ee1,ee2,fft_pad=False)
xcorr2pc = correlate2d(ee2,ee1,fft_pad=True)
xcorr2nc = correlate2d(ee2,ee1,fft_pad=False)

xcorr1pn = correlate2d(ee1,ee2,fft_pad=True,crop=False)
xcorr1nn = correlate2d(ee1,ee2,fft_pad=False,crop=False)
xcorr2pn = correlate2d(ee2,ee1,fft_pad=True,crop=False)
xcorr2nn = correlate2d(ee2,ee1,fft_pad=False,crop=False)

figure(1)
clf()
subplot(221); imshow(abs(acorr1pc)); title("acorr1 pad"); colorbar()
subplot(222); imshow(abs(acorr1nc)); title("acorr1 no pad"); colorbar()
subplot(223); imshow(abs(acorr2pc)); title("acorr2 pad"); colorbar()
subplot(224); imshow(abs(acorr2nc)); title("acorr2 no pad"); colorbar()

figure(4)
clf()
subplot(221); imshow(abs(acorr1pn)); title("acorr1 no crop pad"); colorbar()
subplot(222); imshow(abs(acorr1nn)); title("acorr1 no crop no pad"); colorbar()
subplot(223); imshow(abs(acorr2pn)); title("acorr2 no crop pad"); colorbar()
subplot(224); imshow(abs(acorr2nn)); title("acorr2 no crop no pad"); colorbar()

figure(2)
clf()
subplot(221); imshow(abs(ee1)); title("ee1"); colorbar()
subplot(222); imshow(aconv1p); title("aconv1"); colorbar()
subplot(223); imshow(abs(ee2)); title("ee2"); colorbar()
subplot(224); imshow(aconv2p); title("aconv2"); colorbar()

figure(3)
clf()
subplot(221); imshow(aconv1p); title("aconv1 pad"); colorbar()
subplot(222); imshow(aconv1n); title("aconv1 no pad"); colorbar()
subplot(223); imshow(aconv2p); title("aconv2 pad"); colorbar()
subplot(224); imshow(aconv2n); title("aconv2 no pad"); colorbar()


figure(5)
clf()
subplot(221); imshow(abs(xcorr1pc)); title("xcorr1 pad"); colorbar()
subplot(222); imshow(abs(xcorr1nc)); title("xcorr1 no pad"); colorbar()
subplot(223); imshow(abs(xcorr2pc)); title("xcorr2 pad"); colorbar()
subplot(224); imshow(abs(xcorr2nc)); title("xcorr2 no pad"); colorbar()

figure(6)
clf()
subplot(221); imshow(abs(xcorr1pn)); title("xcorr1 no crop pad"); colorbar()
subplot(222); imshow(abs(xcorr1nn)); title("xcorr1 no crop no pad"); colorbar()
subplot(223); imshow(abs(xcorr2pn)); title("xcorr2 no crop pad"); colorbar()
subplot(224); imshow(abs(xcorr2nn)); title("xcorr2 no crop no pad"); colorbar()

