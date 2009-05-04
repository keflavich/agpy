
procedure pipeline_jhk (filename,outname)
string filename
string outname
bool recrop=yes                 {prompt="Redo cropping / transforming?"}
bool justcrop=no                {prompt="Just crop (no transform?)"}
bool combine=yes                {prompt="Combine images?"}
bool background=yes             {prompt="Run background subtraction?"}
bool calibrate=no               {prompt="Calibrate?"}
bool interactive_background=no  {prompt="Do background subtraction interactively?"}
struct *flist
begin
    struct line
    string filenameJ
    string filenameH
    string filenameK
    real KdeltaL,HdeltaL,JdeltaL,KdeltaP,HdeltaP,JdeltaP
    flist=filename
    delete ( filename+"J" )
    delete ( filename+"H" )
    delete ( filename+"K" )
    print ( "", > filename+"J" )
    print ( "", > filename+"H" )
    print ( "", > filename+"K" )
    if (recrop) {
        while(fscan(flist,line)!=EOF) {
#            imarith (line ,'*', "bmask", line )
            fixpix ( line , "bpm.fits" )
            hedit (line,"DISPAXIS",1,add+,ver-)

            # separate out the individual orders
            imcopy (line+"[*,722:833]", line+"_K.fits", verbose="yes")
            imcopy (line+"[*,533:666]", line+"_H.fits", verbose="yes")
            imcopy (line+"[*,352:537]", line+"_J.fits", verbose="yes")
            filenameJ = line+"_J"
            filenameH = line+"_H"
            filenameK = line+"_K"

            # apply pre-calculated coordinate transforms
            if (justcrop==no) {
                transform ( input=filenameJ , output=filenameJ+"t", fitnames="skylinesJ,starsJ" )
                transform ( input=filenameH , output=filenameH+"t", fitnames="skylinesH,starsH" )
                transform ( input=filenameK , output=filenameK+"t", fitnames="skylinesK,starsK" )
                imcopy (filenameJ+"t[*,53:137]", filenameJ+"tc.fits", verbose="yes")
                imcopy (filenameH+"t[*,18:118]", filenameH+"tc.fits", verbose="yes")
                imcopy (filenameK+"t[*,7:99]", filenameK+"tc.fits", verbose="yes")
                print ( filenameJ+"tc", >> filename+"J" )
                print ( filenameH+"tc", >> filename+"H" )
                print ( filenameK+"tc", >> filename+"K" )
            }
        }
    }

    # combine J/H/K spectra
    if (combine) {
        imcombine ( "@"+filename+"J" , outname+"_J_combine" , combine="median" , scale="mode" )
        imcombine ( "@"+filename+"H" , outname+"_H_combine" , combine="median" , scale="mode" )
        imcombine ( "@"+filename+"K" , outname+"_K_combine" , combine="median" , scale="mode" )
    }


    # calibrate using some sort of standard/sensfunc (I used blackbodies)
    if (calibrate) {
        calibrate ( outname+"_J_combine" , outname+"_J_cal" , sens="sensJ" , ignoreaps+)
        calibrate ( outname+"_H_combine" , outname+"_H_cal" , sens="sensH" , ignoreaps+)
        calibrate ( outname+"_K_combine" , outname+"_K_cal" , sens="sensK" , ignoreaps+)
        imgets ( outname+"_K_cal" , "CD1_1")
        KdeltaL = imgets.value
        imgets ( outname+"_J_cal" , "CD1_1")
        JdeltaL = imgets.value
        imgets ( outname+"_H_cal" , "CD1_1")
        HdeltaL = imgets.value
        imgets ( outname+"_K_cal" , "CD2_2")
        KdeltaP = imgets.value
        imgets ( outname+"_J_cal" , "CD2_2")
        JdeltaP = imgets.value
        imgets ( outname+"_H_cal" , "CD2_2")
        HdeltaP = imgets.value
        magnify ( outname+"_J_cal" , outname+"_J_cal_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
        magnify ( outname+"_H_cal" , outname+"_H_cal_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
        imcombine ( outname+"_K_cal,"+outname+"_H_cal_m,"+outname+"_J_cal_m" , outname+"_JHK_cal" , combine="sum" , offset="wcs" )
        if (background) {
            background ( outname+"_JHK_cal" , outname+"_JHK_cal_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background )
        }
    }
    # background subtract to remove night sky lines
    else {
        if (background) {
            imgets ( outname+"_K_combine" , "CD1_1")
            KdeltaL = imgets.value
            imgets ( outname+"_J_combine" , "CD1_1")
            JdeltaL = imgets.value
            imgets ( outname+"_H_combine" , "CD1_1")
            HdeltaL = imgets.value
            imgets ( outname+"_K_combine" , "CD2_2")
            KdeltaP = imgets.value
            imgets ( outname+"_J_combine" , "CD2_2")
            JdeltaP = imgets.value
            imgets ( outname+"_H_combine" , "CD2_2")
            HdeltaP = imgets.value
            magnify ( outname+"_J_combine" , outname+"_J_combine_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
            magnify ( outname+"_H_combine" , outname+"_H_combine_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
#            background ( outname+"_J_combine" , outname+"_J_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background )
#            background ( outname+"_H_combine" , outname+"_H_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background )
#            background ( outname+"_K_combine" , outname+"_K_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background )
            imcombine ( outname+"_K_combine,"+outname+"_H_combine_m,"+outname+"_J_combine_m" , outname+"_JHK_combine" , combine="sum" , offset="wcs" )
            background ( outname+"_JHK_combine" , outname+"_JHK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background )
        }
    }

end

