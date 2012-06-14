from pylab import *
ion()
import pyfits

def wise_fix_segment(xl,xh,yl,yh,WISE,MSX,minOKratio=None,maxOKratio=None,minOKmsx=None, doplot=True):
    """
    Fix a saturated segment of a WISE or SPITZER image using MSX data
    """

    subMSX = MSX[yl:yh,xl:xh].copy()
    subWISE = WISE[yl:yh,xl:xh].copy()

    mask = np.ones(subMSX.shape, dtype='bool')
    #mask[516:613,73:176] = False
    #mask[610:645,274:346] = False
    #mask[178:585,210:611] = False
    if minOKratio is not None:
        mask[subWISE/subMSX < minOKratio] = False
    if maxOKratio is not None:
        mask[subWISE/subMSX > maxOKratio] = False

    subWISEm = subWISE.copy()
    subWISEm[True-mask] = 0
    mask[np.isnan(subWISEm)] = False
    #subWISEm = subWISE.copy()
    #subWISEm[True-mask] = 0
    #mask[np.isnan(subWISEm)] = False

    if minOKmsx is not None:
        mask[subMSX<minOKmsx] = False

    scalefactor = np.percentile(subWISE[mask]/subMSX[mask],25)
    print "Scalefactor (25%%): %e" % scalefactor
    from agpy import PCA_tools
    m,b = PCA_tools.total_least_squares(subMSX[mask],subWISE[mask])
    print "Scalefactor: %e  offset: %e" % (m,b)


    subWISEm[True-mask] = subMSX[True-mask] * m + b


    newmask = ((MSX * m + b) > WISE * 1.5 )
    if minOKmsx is not None:
        newmask *= (MSX > minOKmsx)
    newmask += np.isnan(WISE)

    import AG_fft_tools
    submask  = newmask[yl:yh,xl:xh]
    subedges = np.sum(np.dstack(np.gradient(submask)),axis=2)
    smsubedges = AG_fft_tools.smooth(subedges)
    smsubWISEm = AG_fft_tools.smooth(subWISEm)
    #edgesmd = subWISEm.copy()
    #edgesmd[True-(smsubedges>1e-2)] = np.nan
    #edgesmd = AG_fft_tools.smooth(edgesmd,interpolate_nan=True)
    subWISEm_fixed = subWISEm.copy()
    subWISEm_fixed[smsubedges>1e-4] = smsubWISEm[smsubedges>1e-4]
    submask += smsubedges > 1e-4

    #WISE[0].data[newmask] = MSX[0].data[newmask]*scalefactor
    newmask[:yl,:] = False
    newmask[yh:,:] = False
    newmask[:,:xl] = False
    newmask[:,xh:] = False
    #WISE[0].data[newmask] = MSX[0].data[newmask]*m+b
    WISE[newmask] = subWISEm_fixed[submask]

    if doplot:
        figure(1)
        clf()
        subplot(211); imshow(subMSX)
        subplot(212); imshow(subWISE)

        figure(2)
        clf()
        H,L,P = hist(subWISE[mask]/subMSX[mask],bins=100)
        vlines(scalefactor,H.min(),H.max(),color='r')
        vlines(m,H.min(),H.max(),color='g')

        figure(3)
        clf()
        loglog(subMSX[mask],subWISE[mask],',')
        xlabel("SubMSX")
        ylabel("SubWISE")
        plot(logspace(-6,-4),logspace(-6,-4)*scalefactor,color='r')
        plot(logspace(-6,-4),logspace(-6,-4)*m+b,color='k')

        figure(4)
        clf()
        imshow(np.log10(subWISEm))

        figure(5)
        clf()
        imshow(newmask)
        title("New Mask")

        figure(8)
        clf()
        imshow(np.isnan(WISE))
        title("Band 3 NANs")

        figure(7)
        clf()
        plot(subWISE[mask]/subMSX[mask],',')
        ylabel('ratio')
        xlabel('WISE/MSX')

        figure(6)
        clf()
        imshow(subWISE/subMSX,vmin=scalefactor/3,vmax=scalefactor*3)
        colorbar()

        figure(9)
        clf()
        imshow(np.log10(WISE[yl:yh,xl:xh]))

    return WISE



