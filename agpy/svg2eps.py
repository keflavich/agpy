#!/usr/bin/env python
"""
Convert an SVG to an EPS file 
http://stackoverflow.com/questions/9793120/how-to-covert-svg-to-eps-in-ghostscript
"""
import optparse
import os
import re
import tempfile


if __name__ == "__main__":

    parser=optparse.OptionParser()
    parser.add_option("--remove","-r",help="Remove .ps file after conversion?",default=False,action="store_true")
#    parser.add_option("--multipage","-m",help="Convert to multiple png files?",default=False,action="store_true")
    parser.add_option("--silent","-s",help="Be quiet? Default False",default=False,action="store_true")
    parser.add_option("--verbose","-v",help="Be loud? Default True",default=1)
    parser.add_option("--silence_gs",help="Silence Ghostscript?  Default True",default=True)
#    parser.add_option("--resolution",help="Resolution (pixels per inch) of png.  Default 300",default=300)
    parser.add_option("--noepscrop",help="No EPS crop?  Default False",default=False,action='store_true')
    parser.add_option("--outfile",help="Outfile name?",default=None)

    options,args = parser.parse_args()

    verbose = not(options.silent) and options.verbose

    if verbose > 1:
        print "Args: ",args
        print "Options: ",options

    for filename in args:

        if options.noepscrop: epscrop=""
        else: epscrop = "-dEPSCrop"


#        if options.multipage:
#            outfile = re.sub("\.e?ps","_%d.png",filename)
#        else:
#            outfile = re.sub("\.e?ps",".png",filename)
        if options.outfile is None:
            outfile = re.sub("\.svg$",".eps",filename)
        else:
            outfile = options.outfile
        ps_outfile = re.sub("\.eps",".ps",outfile)

        print "outfile: %s  ps_outfile: %s" % (outfile,ps_outfile)

        command1 = "gsvg -dNOPAUSE -sDEVICE=ps2write -sOutputFile=%s" % (ps_outfile)
        command2 = "ps2eps -f %s" % (ps_outfile)
        #command2 = "gs   -dBATCH -dNOPAUSE -sDEVICE=epswrite %s" % (epscrop)
#        command = "gs -dBATCH -sDEVICE=png16m -r%i %s -dNOPAUSE" % (options.resolution,epscrop)

        command1 += " %s" % filename
        #command2 += " -sOutputFile=%s %s" % (outfile, ps_outfile)

        if options.silence_gs:
            command1 += " > /dev/null"
            command2 += " > /dev/null"

        if verbose: 
            print command1
            print command2
        status = os.system(command1)
        status = os.system(command2)

        if options.remove:
            if verbose > 1: print "rm %s" % filename
            os.remove(filename)





#gsvg \
#    -dNOPAUSE \
#    -sDEVICE=ps2write \
#    -sOutputFile=my.ps \
#    my.svg 
#
#gs \
#    -dNOPAUSE \
#    -sDEVICE=epswrite \
#    -sOutputFile=my.eps \
#    my.ps 
