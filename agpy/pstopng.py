import optparse
import os
import re


if __name__ == "__main__":

    parser=optparse.OptionParser()
    parser.add_option("--remove","-r",help="Remove .ps file after conversion?",default=False,action="store_true")
    parser.add_option("--multipage","-m",help="Convert to multiple png files?",default=False,action="store_true")
    parser.add_option("--silent","-s",help="Be quiet? Default False",default=False,action="store_true")
    parser.add_option("--verbose","-v",help="Be loud? Default True",default=1)
    parser.add_option("--silence_gs",help="Silence Ghostscript?  Default True",default=True)
    parser.add_option("--resolution",help="Resolution (pixels per inch) of png.  Default 300",default=300)
    parser.add_option("--noepscrop",help="No EPS crop?  Default False",default=False,action='store_true')

    options,args = parser.parse_args()

    verbose = not(options.silent) and options.verbose

    if verbose > 1:
        print "Args: ",args
        print "Options: ",options

    for filename in args:

        if options.noepscrop: epscrop=""
        else: epscrop = "-dEPSCrop"

        command = "gs -dBATCH -sDEVICE=png16m -r%i %s -dNOPAUSE" % (options.resolution,epscrop)

        if options.multipage:
            outfile = re.sub("\.e?ps","_%d.png",filename)
        else:
            outfile = re.sub("\.e?ps",".png",filename)

        command += " -sOutputFile=%s %s" % (outfile,filename)

        if options.silence_gs:
            command += " > /dev/null"

        if verbose: print command
        status = os.system(command)

        if options.remove:
            if verbose > 1: print "rm %s" % filename
            os.remove(filename)


