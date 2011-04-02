#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
import region_photometry

def region_photometry_files(regfile,filelist):
    
    allprops = []
    for fn in filelist:
        allprops.append([fn] + region_photometry.region_photometry(regfile,fn,doprint=False))

    return allprops


if __name__ == "__main__":

    import sys
    regionfilename = sys.argv[1]
    fitsfiles = sys.argv[2:]
    print "region: ",regionfilename
    print "fitsfiles: ",fitsfiles

    props = region_photometry_files(regionfilename,fitsfiles)
    print props

    sys.exit()
