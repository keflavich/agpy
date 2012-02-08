
procedure pipeline_jhk (filename,outname)
string filename
string outname
string procdir='./'             {prompt="Output directory for intermediate files"}
string rawdir='./'              {prompt="Directory in which the raw files are stored"}
string dark=''                  {prompt="Filename of dark (assumes current dir)"}
bool darksub=yes                {prompt="Do dark subtraction? (must specify dark!)"}
bool recrop=yes                 {prompt="Redo cropping / transforming?"}
bool justcrop=no                {prompt="Just crop (no transform?)"}
bool combine=yes                {prompt="Combine images?"}
bool background=yes             {prompt="Run background subtraction?"}
bool calibrate=no               {prompt="Calibrate?"}
bool interactive_background=no  {prompt="Do background subtraction interactively?"}
bool calibrator=no              {prompt="Is the source a calibrator?"}
bool clobber=no                 {prompt="Overwrite existing files?"}
bool caleach=no                 {prompt="Calibrate EACH file? (requires pre-existing calibrator data)"}
string sample="*"               {prompt="Selection region for background task"}
real magJ=0.0                   {prompt="J Magnitude of calibrator star"}
real magH=0.0                   {prompt="H Magnitude of calibrator star"}
real magK=0.0                   {prompt="K Magnitude of calibrator star"}
real teff=10000                 {prompt="Effective temperature of calibrator (required for calibrator=yes)"}
struct *flist
begin
    struct line
    string filenameJ
    string filenameH
    string filenameK
    string ds
    real KdeltaL,HdeltaL,JdeltaL,KdeltaP,HdeltaP,JdeltaP
    flist=filename
    delete ( filename+"J" )
    delete ( filename+"H" )
    delete ( filename+"K" )
    print ( "", > filename+"J" )
    print ( "", > filename+"H" )
    print ( "", > filename+"K" )

    onedspec
    twodspec
    apextract
    longslit

    if (recrop) {
        while(fscan(flist,line)!=EOF) {
#            imarith (line ,'*', "bmask", line )
            fixpix ( rawdir+line , "bpm.fits" )
            hedit ( rawdir+line,"DISPAXIS",1,add+,ver-)
            if (darksub) { 
                imarith ( rawdir+line, '-', dark, procdir+line+'_darksub.fits' ) 
                ds = '_darksub' 
                print ( "Subtracted dark "+dark+" from "+line )
                # separate out the individual orders
                imcopy (procdir+line+"_darksub[*,722:833]", procdir+line+"_K.fits", verbose="yes")
                imcopy (procdir+line+"_darksub[*,533:666]", procdir+line+"_H.fits", verbose="yes")
                imcopy (procdir+line+"_darksub[*,352:537]", procdir+line+"_J.fits", verbose="yes")
            } else {
                imcopy (rawdir+line+"[*,722:833]", procdir+line+"_K.fits", verbose="yes")
                imcopy (rawdir+line+"[*,533:666]", procdir+line+"_H.fits", verbose="yes")
                imcopy (rawdir+line+"[*,352:537]", procdir+line+"_J.fits", verbose="yes")
            }

            filenameJ = line+"_J"
            filenameH = line+"_H"
            filenameK = line+"_K"

            # apply pre-calculated coordinate transforms
            if (justcrop==no) {
                transform ( input=procdir+filenameJ , output=procdir+filenameJ+"t", fitnames="skylinesJ,starsJ" )
                transform ( input=procdir+filenameH , output=procdir+filenameH+"t", fitnames="skylinesH,starsH" )
                transform ( input=procdir+filenameK , output=procdir+filenameK+"t", fitnames="skylinesK,starsK" )
                imcopy (procdir+filenameJ+"t[*,53:137]", procdir+filenameJ+"tc.fits", verbose="yes" )
                imcopy (procdir+filenameH+"t[*,18:118]", procdir+filenameH+"tc.fits", verbose="yes" )
                imcopy (procdir+filenameK+"t[*,7:99]",   procdir+filenameK+"tc.fits", verbose="yes" )
                print ( procdir+filenameJ+"tc", >> filename+"J" )
                print ( procdir+filenameH+"tc", >> filename+"H" )
                print ( procdir+filenameK+"tc", >> filename+"K" )
            }

        }
    }

    # calibrate using some sort of standard/sensfunc (I used blackbodies)
    if (caleach) {
        while(fscan(flist,line)!=EOF) {
            filenameJ = line+"_J"
            filenameH = line+"_H"
            filenameK = line+"_K"

            if (clobber) {
                delete ( procdir+filenameJ+"tc_cal.fits" )
                delete ( procdir+filenameH+"tc_cal.fits" )
                delete ( procdir+filenameK+"tc_cal.fits" )
            }

            calibrate ( procdir+filenameJ+"tc" , procdir+filenameJ+"tc_cal" , sens="sensJ" , ignoreaps+)
            calibrate ( procdir+filenameH+"tc" , procdir+filenameH+"tc_cal" , sens="sensH" , ignoreaps+)
            calibrate ( procdir+filenameK+"tc" , procdir+filenameK+"tc_cal" , sens="sensK" , ignoreaps+)
            imgets ( procdir+filenameK+"tc_cal" , "CD1_1")
            KdeltaL = imgets.value
            imgets ( procdir+filenameJ+"tc_cal" , "CD1_1")
            JdeltaL = imgets.value
            imgets ( procdir+filenameH+"tc_cal" , "CD1_1")
            HdeltaL = imgets.value
            imgets ( procdir+filenameK+"tc_cal" , "CD2_2")
            KdeltaP = imgets.value
            imgets ( procdir+filenameJ+"tc_cal" , "CD2_2")
            JdeltaP = imgets.value
            imgets ( procdir+filenameH+"tc_cal" , "CD2_2")
            HdeltaP = imgets.value

            if (clobber) {
                delete ( procdir+filenameJ+"tc_cal_m.fits" )
                delete ( procdir+filenameH+"tc_cal_m.fits" )
                delete ( procdir+filename+"_JHK_cal.fits"  )
            }

            magnify ( procdir+filenameJ+"tc_cal" , procdir+filenameJ+"tc_cal_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
            magnify ( procdir+filenameH+"tc_cal" , procdir+filenameH+"tc_cal_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
            imcombine ( procdir+filenameK+"tc_cal,"+procdir+filenameH+"tc_cal_m,"+procdir+filenameJ+"tc_cal_m" , procdir+line+"_JHK_cal" , combine="sum" , offset="wcs" )
            if (background) {
                if (clobber) {
                    delete ( procdir+line+"_JHK_cal_backsub.fits" )
                }
                background ( procdir+line+"_JHK_cal" , procdir+line+"_JHK_cal_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
            }
        }
    }

    # combine J/H/K spectra
    if (combine) {
        imcombine ( "@"+filename+"J" , procdir+outname+"_J_combine" , combine="median" , scale="mode" )
        imcombine ( "@"+filename+"H" , procdir+outname+"_H_combine" , combine="median" , scale="mode" )
        imcombine ( "@"+filename+"K" , procdir+outname+"_K_combine" , combine="median" , scale="mode" )
    }


    # calibrate using some sort of standard/sensfunc (I used blackbodies)
    if (calibrate) {
        calibrate ( procdir+outname+"_J_combine" , procdir+outname+"_J_cal" , sens="sensJ" , ignoreaps+)
        calibrate ( procdir+outname+"_H_combine" , procdir+outname+"_H_cal" , sens="sensH" , ignoreaps+)
        calibrate ( procdir+outname+"_K_combine" , procdir+outname+"_K_cal" , sens="sensK" , ignoreaps+)
        imgets ( procdir+outname+"_K_cal" , "CD1_1")
        KdeltaL = imgets.value
        imgets ( procdir+outname+"_J_cal" , "CD1_1")
        JdeltaL = imgets.value
        imgets ( procdir+outname+"_H_cal" , "CD1_1")
        HdeltaL = imgets.value
        imgets ( procdir+outname+"_K_cal" , "CD2_2")
        KdeltaP = imgets.value
        imgets ( procdir+outname+"_J_cal" , "CD2_2")
        JdeltaP = imgets.value
        imgets ( procdir+outname+"_H_cal" , "CD2_2")
        HdeltaP = imgets.value
        magnify ( procdir+outname+"_J_cal" , procdir+outname+"_J_cal_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
        magnify ( procdir+outname+"_H_cal" , procdir+outname+"_H_cal_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
        imcombine ( procdir+outname+"_K_cal,"+procdir+outname+"_H_cal_m,"+procdir+outname+"_J_cal_m" , procdir+outname+"_JHK_cal" , combine="sum" , offset="wcs" )
        if (background) {
            background ( procdir+outname+"_JHK_cal" , procdir+outname+"_JHK_cal_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
        }
    }
    # background subtract to remove night sky lines
    else {
        if (background) {
            imgets ( procdir+outname+"_K_combine" , "CD1_1")
            KdeltaL = imgets.value
            imgets ( procdir+outname+"_J_combine" , "CD1_1")
            JdeltaL = imgets.value
            imgets ( procdir+outname+"_H_combine" , "CD1_1")
            HdeltaL = imgets.value
            imgets ( procdir+outname+"_K_combine" , "CD2_2")
            KdeltaP = imgets.value
            imgets ( procdir+outname+"_J_combine" , "CD2_2")
            JdeltaP = imgets.value
            imgets ( procdir+outname+"_H_combine" , "CD2_2")
            HdeltaP = imgets.value
            magnify ( procdir+outname+"_J_combine" , procdir+outname+"_J_combine_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
            magnify ( procdir+outname+"_H_combine" , procdir+outname+"_H_combine_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
            if (calibrator) {
                background ( procdir+outname+"_J_combine" , procdir+outname+"_J_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                background ( procdir+outname+"_H_combine" , procdir+outname+"_H_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                background ( procdir+outname+"_K_combine" , procdir+outname+"_K_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                apall ( procdir+outname+"_J_backsub" , interactive="no" )
                apall ( procdir+outname+"_H_backsub" , interactive="no" )
                apall ( procdir+outname+"_K_backsub" , interactive="no" )
                standard ( procdir+outname+"_J_backsub.ms" , output="stdJ" , star_name="J" , caldir="onedstds$blackbody/" , mag=magJ, teff=teff , magband="J" )
                standard ( procdir+outname+"_H_backsub.ms" , output="stdH" , star_name="H" , caldir="onedstds$blackbody/" , mag=magH, teff=teff , magband="H" )
                standard ( procdir+outname+"_K_backsub.ms" , output="stdK" , star_name="K" , caldir="onedstds$blackbody/" , mag=magK, teff=teff , magband="K" )
                sensfunc ( "stdJ" , "sensJ" , answer="YES" )
                sensfunc ( "stdH" , "sensH" , answer="YES" )
                sensfunc ( "stdK" , "sensK" , answer="YES" )
            }
            imcombine ( procdir+outname+"_K_combine,"+procdir+outname+"_H_combine_m,"+procdir+outname+"_J_combine_m" , procdir+outname+"_JHK_combine" , combine="sum" , offset="wcs" )
            background ( procdir+outname+"_JHK_combine" , procdir+outname+"_JHK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
        }
    }

end

