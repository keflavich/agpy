#!/bin/bash

# Wrapper to mosaic a *subset* of images in the current directory
#
# Usage:
# mGetHdr template_fitsfile.fit mosaic.hdr
# montage outfile=l089_reference_montage.fits 0*_indiv13pca*_map01.fits combine=median &
#
# Keyword parameters:
#   combine - 'median', 'average', or 'sum'
#   header  - Filename of fits file from which to create mosaic.hdr
#   outfile - Output fits filename
#

#export PATH=$PATH:/usr/local/bin/montage

origdir=`pwd`

debug=0
copy="F"

if [ $# -eq 0 ]
then
    echo montage wrapper
    echo Usage:
    echo mGetHdr template_fitsfile.fit mosaic.hdr
    echo montage outfile=l089_reference_montage.fits 0\*_indiv13pca\*_map01.fits combine='median' \&
    echo
    echo Keyword parameters:
    echo "  combine - 'median', 'average', or 'sum'  (note that this has to be \"combine='median'\")"
    echo "  exact_size - True or False"
    echo "  header  - Filename of fits file from which to create mosaic.hdr (or the .hdr file)"
    echo "             (if header is not specified, assumed to be mosaic.hdr in pwd)"
    echo "  outfile - Output fits filename"
    echo "  copy    - copy files instead of hardlinking"
    exit
else
    for ii in $*
    do
        if [ ${ii%=*} == 'header' ]
        then
            hdrfile=${ii#*=}
            if [ ${hdrfile#*.} == 'hdr' ]
            then
              headerfile=$hdrfile
            elif [ ${hdrfile#*.} == 'fits' ]
            then
              /usr/local/bin/montage/mGetHdr $hdrfile mosaic.hdr
              headerfile=mosaic.hdr
            else
              echo "ERROR: $ii is not a valid header setting"
              headerfile="does_not_exist_I_hope.noooo"
            fi
        elif [ ${ii} == '-copy' ]
        then
            copy="T"
        elif [ ${ii%=*} == 'outfile' ]
        then
            outfile=${ii#*=}
        elif [ `echo $ii | grep =` ] 
        then
            params="$params,${ii%=*}=${ii#*=}"
        # note that you can say "header=blah.fits" because it is checked for earlier
        elif [ $ii == "debug" ] 
        then
            debug=1
        elif [ `echo $ii | grep ".fits"` ] 
        then
            files=( ${files[@]} $ii )
        fi
    done
fi
echo FILES: ${files[@]} 
echo NFILES: ${#files} 
echo Extra Parameters: $params
if [ ${#files} -gt 0 ] 
then
    echo "Creating temporary directory and hard-linking (not sym-linking) all files into it"
    #echo "Creating temporary directory and sym-linking all files into it"
    mkdir tmp
    if [ $copy == "T" ] 
    then 
        cp ${files[@]} tmp/
        cp $headerfile tmp/
    else
        #echo ln ${files[@]} tmp/
        ln ${files[@]} tmp/
        #echo ln $headerfile tmp/
        ln $headerfile tmp/
    fi
    cd tmp/
    #ls -lh
fi

if [ -f $headerfile ] 
then 
    echo "$headerfile exists, continuing"

    if [ $debug == 1 ] 
    then
      cmd="/usr/bin/env ipython -pdb -c "
    else
      cmd="/usr/bin/env python -c "
    fi

    montagecmd="import tempfile; tempfile.tempdir='/Volumes/disk4/var/tmp'; import montage; montage.wrappers.mosaic('$dir','$dir/mosaic',header='$dir/$headerfile'$params)"
    dir=`pwd`
    echo $cmd $montagecmd
    $cmd $montagecmd 
    if [ $? != 0 ]
    then
      echo "There was a python error!"
    fi

    cd $origdir
    if [ -d tmp ]
    then
        if [ $outfile ]
        then
            mv tmp/mosaic/mosaic.fits $outfile
        else
            mv tmp/mosaic mosaic
        fi
        rm -r tmp
    fi

else
    echo "$headerfile does not exist.  Quitting."
    cd $origdir
fi

exit
