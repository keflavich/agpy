
procedure pipeline_dis ()
string filename
string outname
string letter
string dispersion,centerwavelength,platescale,angle,exptime
string imtype,detector,objectname,lamp,obj
struct *flist
begin
    int red,blue,listind,wherespace
#    mkdir ( "pipeline")
    struct line
    if (access("pipeline/allfits")==no) {
        files ( "*.fits" , > "pipeline/allfits" )
        flist="pipeline/allfits"
        delete ( "pipeline/bluecallist" )
        delete ( "pipeline/blueobjlist" )
        delete ( "pipeline/bluebiaslist" )
        while(fscan(flist,line)!=EOF) {
            print ( line ) 
            imgets(line,"DISPDW")
            dispersion = imgets.value
            imgets(line,"DISPWC")
            centerwavelength = imgets.value
            imgets(line,"PIXSCAL2")
            platescale = imgets.value
            imgets(line,"IMAGETYP")
            imtype = imgets.value
            imgets(line,"DETECTOR")
            detector = imgets.value
            imgets(line,"OBJNAME")
            objectname = imgets.value
            imgets(line,"EXPTIME")
            exptime = imgets.value
            imgets(line,"OBJANGLE")
            angle = imgets.value
            if (angle!=0) { print(angle) | scanf("%5f",angle) } 
            imgets(line,"LAMP")
            lamp = imgets.value

            wherespace=stridx(" ",objectname)
            while (wherespace > 0) {
                objectname=substr(objectname,1,wherespace-1)+substr(objectname,wherespace+1,strlen(objectname))
                wherespace=stridx(" ",objectname)
            }
        
            if (detector=="red") { 
                print (dispersion) | scanf("%f",dispersion)
                dispersion=dispersion*-1.0 }

            print ( "Center Wavelength: " , centerwavelength , " dispersion " , dispersion )
            hedit(line, "CD1_1",  0, add=yes , verify=no , addonly=no)
            hedit(line, "CRVAL1", 0, add=yes , verify=no , addonly=no)
            hedit(line, "CD1_1", dispersion,   verify=no  )
            hedit(line, "CRVAL1", centerwavelength,  verify=no )

            hedit(images=line, fields="CRPIX1", value=1024 , add+ , verify- )
            hedit(images=line, fields="CRPIX2", value=512 , add+ , verify- )
            hedit(images=line, fields="CRVAL2", value=512 , add+ , verify- )
            hedit(images=line, fields="CD2_2", value=1 , add+ , verify- )
            hedit(images=line, fields="CD1_2", value=0 , add+ , verify- )
            hedit(images=line, fields="CD2_1", value=0 , add+ , verify- )

            if (lamp == "He") {
                if (detector=="blue") {print ( line , >> "pipeline/blueHelist" )}
                if (detector=="red")  {print ( line , >> "pipeline/redHelist"  )}}
            else if (lamp == "Ar") {
                if (detector=="blue") {print ( line , >> "pipeline/blueArlist" )}
                if (detector=="red")  {print ( line , >> "pipeline/redArlist"  )}}
            else if (lamp == "Ne") {
                if (detector=="blue") {print ( line , >> "pipeline/blueNelist" )}
                if (detector=="red")  {print ( line , >> "pipeline/redNelist"  )}}
            else if (imtype == "zero") {
                if (detector=="blue") {print ( line , >> "pipeline/blueBias" )}
                if (detector=="red")  {print ( line , >> "pipeline/redBias"  )}}
            else {
                if (detector=="blue") {print ( line , >> "pipeline/blue"+objectname+angle+"list" )}
                if (detector=="red")  {print ( line , >> "pipeline/red"+objectname+angle+"list"  )}}
        }
    }

    if (access("pipeline/listoflists")==yes) { delete("pipeline/listoflists") }
    files ( "pipeline/*list" , > "pipeline/listoflists" )
    flist="pipeline/listoflists"

    imcombine ( "@pipeline/BlueBias" , "pipeline/BlueBias.fits" , combine="median", scale="mode" )
    imcombine ( "@pipeline/RedBias" , "pipeline/RedBias.fits" , combine="median", scale="mode" )

    while(fscan(flist,line)!=EOF) {
        red = strstr("red",line)
        blue = strstr("blue",line)
        listind = strstr("list",line)
#        sort (line) | unique > $line
        if (red > 1) {
            obj = substr(line,red+3,listind-1)
            print ( line,red,blue,listind,obj )
#            imarith ( "@"+line  , "-", "pipeline/RedBias.fits", "@"+line)
            imcombine ( "@"+line,  "pipeline/Red"+obj+".fits" , combine="median", scale="mode" )
            print ("Red"+obj,>>"pipeline/RedObjectList")
        } else if (blue > 1) {
            obj = substr(line,blue+4,listind-1)
            print ( line,red,blue,listind,obj )
#            imarith ( "@"+line  , "-", "pipeline/BlueBias.fits", "@"+line)
            imcombine ( "@"+line , "pipeline/Blue"+obj+".fits" , combine="median", scale="mode" )
            print ("Blue"+obj,>>"pipeline/BlueObjectList")
        }

    }
    
    chdir ( "pipeline" )
    files ("RedHe.fits" , >> "RedHeNeAr")
    files ("RedNe.fits" , >> "RedHeNeAr")
    files ("RedAr.fits" , >> "RedHeNeAr")
    files ("BlueHe.fits",  >> "BlueHeNeAr")
    files ("BlueNe.fits",  >> "BlueHeNeAr")
    files ("BlueAr.fits",  >> "BlueHeNeAr")
    imcombine ( "@BlueHeNeAr" , "BlueHeNeAr.fits" , combine="sum")
    imcombine ( "@RedHeNeAr" , "RedHeNeAr.fits" , combine="sum" )

    imgets("RedHeNeAr","CD1_1")
    dispersion = imgets.value
    imgets("RedHeNeAr","CRVAL1")
    centerwavelength = imgets.value
    print ("Red auto with center: ",centerwavelength," and dispersion ",dispersion)
    autoidentify("RedHeNeAr", centerwavelength, dispersion, coordlist="linelists$idhenear.dat")
    reidentify("RedHeNeAr","RedHeNeAr", answer="NO", interactive-, nlost=5, verbose+)
    fitcoords("RedHeNeAr",xorder=4,yorder=3,interactive-)

    imgets("BlueHeNeAr","CD1_1")
    dispersion = imgets.value
    imgets("BlueHeNeAr","CRVAL1")
    centerwavelength = imgets.value
    autoidentify("BlueHeNeAr", centerwavelength, dispersion, coordlist="linelists$henearhres.dat")
    reidentify("BlueHeNeAr","BlueHeNeAr", answer="NO", interactive-, nlost=20, verbose+)
    fitcoords("BlueHeNeAr",xorder=4,yorder=3,interactive-)

    flist="RedObjectList"
    while(fscan(flist,line)!=EOF) {
        transform(line,line+"_trans",fitnames="RedHeNeAr")
        background(line+"_trans[1:2048,206:958]",line+"_backsub",axis=2,inter-,niter=3,high=2,low=2,order=3,sample="1:2048,206:958")
    }
    flist="BlueObjectList"
    while(fscan(flist,line)!=EOF) {
        transform(line,line+"_trans",fitnames="BlueHeNeAr")
        background(line+"_trans[1:2048,98:802]",line+"_backsub",axis=2,inter-,niter=3,high=2,low=2,order=3,sample="1:2048,98:802")
    }

    chdir("..")
end
    
