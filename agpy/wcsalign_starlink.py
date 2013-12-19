#!/bin/env python
from __future__ import print_function

#  Re-projects a 2D FITS image so that it is aligned pixel-for-pixel with
#  another reference FITS image.

#  Usage:
#    wcsalign.py <infile> <reffile> <outfile>
#
#    <infile> = Input fits file
#    <reffile> = Reference fits file
#    <outfile> = Output fits file

#  Notes:
#     - Requires pyast version 2.3

try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
import sys

def main():
    #  Open the input and reference FITS files using fits. A list of the HDUs in each FITS
    #  file is returned.
    hdu_list_in = pyfits.open(sys.argv[1])
    hdu_list_ref = pyfits.open(sys.argv[2])

    wcsalign(hdu_list_in, hdu_list_ref)

def wcsalign(hdu_in, header, outname=None, clobber=False):
    """
    Align one FITS image to a specified header

    Requires pyast.

    Parameters
    ----------
    hdu_in: fits.PrimaryHDU
        The HDU to reproject (must have header & data)
    header: fits.Header
        The target header to project to
    outname: str (optional)
        The filename to write to.
    clobber: bool
        Overwrite the file `outname` if it exists

    Returns
    -------
    The reprojected fits.PrimaryHDU

    Credits
    -------
    Written by David Berry and adapted to functional form by Adam Ginsburg
    (adam.g.ginsburg@gmail.com)
    """
    import starlink.Ast as Ast
    import starlink.Atl as Atl
  
    #  Create objects that will transfer FITS header cards between an AST
    #  FitsChan and the fits header describing the primary HDU of the
    #  supplied FITS file.
    adapter_in = Atl.PyFITSAdapter(hdu_in)
    hdu_ref = pyfits.PrimaryHDU(header=header)
    adapter_ref = Atl.PyFITSAdapter(hdu_ref)
     
    #  Create a FitsChan for each and use the above adapters to copy all the header
    #  cards into it.
    fitschan_in = Ast.FitsChan(adapter_in, adapter_in)
    fitschan_ref = Ast.FitsChan(adapter_ref, adapter_ref)
     
    #  Get the flavour of FITS-WCS used by the header cards currently in the
    #  input FITS file. This is so that we can use the same flavour when we write
    #  out the modified WCS.
    encoding = fitschan_in.Encoding
     
    #  Read WCS information from the two FitsChans. Additionally, this removes
    #  all WCS information from each FitsChan. The returned wcsinfo object
    #  is an AST FrameSet, in which the current Frame describes WCS coordinates
    #  and the base Frame describes pixel coodineates. The FrameSet includes a
    #  Mapping that specifies the transformation between the two Frames.
    wcsinfo_in = fitschan_in.read()
    wcsinfo_ref = fitschan_ref.read()
     
    #  Check that the input FITS header contained WCS in a form that can be
    #  understood by AST.
    if wcsinfo_in is None:
        raise ValueError("Failed to read WCS information from {0}".format(hdu_in))
     
    #  This is restricted to 2D arrays, so check theinput  FITS file has 2 pixel
    #  axes (given by Nin) and 2 WCS axes (given by Nout).
    elif wcsinfo_in.Nin != 3 or wcsinfo_in.Nout != 3:
        raise ValueError("{0} is not 3-dimensional".format(hdu_in))
     
    #  Check the reference FITS file in the same way.
    elif wcsinfo_ref is None:
        raise ValueError("Failed to read WCS information from {0}".format(hdu_ref))
     
    elif wcsinfo_ref.Nin != 3 or wcsinfo_ref.Nout != 3:
        raise ValueError("{0} is not 2-dimensional".format(hdu_ref))
     
    #  Proceed if the WCS information was read OK.
 
    #  Attempt to get a mapping from pixel coords in the input FITS file to pixel
    #  coords in the reference fits file, with alignment occuring by preference in
    #  the current WCS frame. Since the pixel coordinate frame will be the base frame
    #  in each Frameset, we first invert the FrameSets. This is because the Convert method
    #  aligns current Frames, not base frames.
    wcsinfo_in.invert()
    wcsinfo_ref.invert()
    alignment_fs = wcsinfo_in.convert(wcsinfo_ref)

    #  Invert them again to put them back to their original state (i.e.
    #  base frame = pixel coords, and current Frame = WCS coords).
    wcsinfo_in.invert()
    wcsinfo_ref.invert()

    #  Check alignment was possible.
    if alignment_fs is None:
        raise Exception("Cannot find a common coordinate system shared by {0} and {1}".format(hdu_in,hdu_ref))
      
    else:
        #  Get the lower and upper bounds of the input image in pixel indices.
        #  All FITS arrays by definition have lower pixel bounds of [1,1] (unlike
        #  NDFs). Note, unlike fits AST uses FITS ordering for storing pixel axis
        #  values in an array (i.e. NAXIS1 first, NAXIS2 second, etc).
        lbnd_in = [1, 1, 1]
        ubnd_in = [fitschan_in["NAXIS1"], fitschan_in["NAXIS2"], fitschan_in["NAXIS3"]]
 
        #  Find the pixel bounds of the input image within the pixel coordinate
        #  system of the reference fits file.
        (lb1, ub1, xl, xu) = alignment_fs.mapbox(lbnd_in, ubnd_in, 1)
        (lb2, ub2, xl, xu) = alignment_fs.mapbox(lbnd_in, ubnd_in, 2)
 
        #  Calculate the bounds of the output image.
        lbnd_out = [int(lb1), int(lb2)]
        ubnd_out = [int(ub1), int(ub2)]
      
        #  Unlike NDFs, FITS images cannot have an arbitrary pixel origin so
        #  we need to ensure that the bottom left corner of the input image
        #  gets mapped to pixel [1,1] in the output. To do this we, extract the
        #  mapping from the alignment FrameSet and add on a ShiftMap (a mapping
        #  that just applies a shift to each axis).
        shift = [1 - lbnd_out[0],
                 1 - lbnd_out[1]]
      
        alignment_mapping = alignment_fs.getmapping()
        shiftmap = Ast.ShiftMap(shift)
        total_map = Ast.CmpMap(alignment_mapping, shiftmap)
      
        #  Modify the pixel bounds of the output image to take account of this
        #  shift of origin.
        lbnd_out[0] += shift[0]
        lbnd_out[1] += shift[1]
        ubnd_out[0] += shift[0]
        ubnd_out[1] += shift[1]
 
        #  Get the value used to represent missing pixel values
        if "BLANK" in fitschan_in:
            badval = fitschan_in["BLANK"]
            flags = Ast.USEBAD
        else:
            badval = 0
            flags = 0
 
        # Resample the data array using the above mapping.
        # total_map was pixmap; is this right?
        (npix, out, out_var) = total_map.resample(lbnd_in, ubnd_in,
                                                  hdu_in.data, None,
                                                  Ast.LINEAR, None, flags,
                                                  0.05, 1000, badval, lbnd_out,
                                                  ubnd_out, lbnd_out, ubnd_out)
 
        #  Store the aligned data in the primary HDU, and update the NAXISi keywords
        #  to hold the number of pixels along each edge of the rotated image.
        hdu_in.data = out
        fitschan_in["NAXIS1"] = ubnd_out[0] - lbnd_out[0] + 1
        fitschan_in["NAXIS2"] = ubnd_out[1] - lbnd_out[1] + 1

    #  The WCS to store in the output is the same as the reference WCS
    #  except for the extra shift of origin. So use the above shiftmap to
    #  remap the pixel coordinate frame in the reference WCS FrameSet. We
    #  can then use this FrameSet as the output FrameSet.
    wcsinfo_ref.remapframe(Ast.BASE, shiftmap)
 
    #  Attempt to write the modified WCS information to the primary HDU (i.e.
    #  convert the FrameSet to a set of FITS header cards stored in the
    #  FITS file). Indicate that we want to use original flavour of FITS-WCS.
    fitschan_in.Encoding = encoding
    fitschan_in.clear('Card')

    if fitschan_in.write(wcsinfo_ref) == 0:
        raise Exception("Failed to convert the aligned WCS to Fits-WCS")
 
    #  If successful, force the FitsChan to copy its contents into the
    #  fits header, then write the changed data and header to the output
    #  FITS file.
    else:
        fitschan_in.writefits()

    if outname is not None:
        hdu_in.writeto(outname, clobber=clobber)
    
    return hdu_in


def wcsalign_cube(hdu_in, header, outname=None, clobber=False):
    """
    Align one FITS image *cube* to a specified header

    Requires pyast.

    Parameters
    ----------
    hdu_in: fits.PrimaryHDU
        The HDU to reproject (must have header & data)
    header: fits.Header
        The target header to project to
    outname: str (optional)
        The filename to write to.
    clobber: bool
        Overwrite the file `outname` if it exists

    Returns
    -------
    The reprojected fits.PrimaryHDU

    Credits
    -------
    Written by David Berry and adapted to functional form by Adam Ginsburg
    (adam.g.ginsburg@gmail.com)
    """
    import starlink.Ast as Ast
    import starlink.Atl as Atl
  
    #  Create objects that will transfer FITS header cards between an AST
    #  FitsChan and the fits header describing the primary HDU of the
    #  supplied FITS file.
    adapter_in = Atl.PyFITSAdapter(hdu_in)
    hdu_ref = pyfits.PrimaryHDU(header=header)
    adapter_ref = Atl.PyFITSAdapter(hdu_ref)
     
    #  Create a FitsChan for each and use the above adapters to copy all the header
    #  cards into it.
    fitschan_in = Ast.FitsChan(adapter_in, adapter_in)
    fitschan_ref = Ast.FitsChan(adapter_ref, adapter_ref)
     
    #  Get the flavour of FITS-WCS used by the header cards currently in the
    #  input FITS file. This is so that we can use the same flavour when we write
    #  out the modified WCS.
    encoding = fitschan_in.Encoding
     
    #  Read WCS information from the two FitsChans. Additionally, this removes
    #  all WCS information from each FitsChan. The returned wcsinfo object
    #  is an AST FrameSet, in which the current Frame describes WCS coordinates
    #  and the base Frame describes pixel coodineates. The FrameSet includes a
    #  Mapping that specifies the transformation between the two Frames.
    wcsinfo_in = fitschan_in.read()
    wcsinfo_ref = fitschan_ref.read()
     
    #  Check that the input FITS header contained WCS in a form that can be
    #  understood by AST.
    if wcsinfo_in is None:
        raise ValueError("Failed to read WCS information from {0}".format(hdu_in))
     
    #  This is restricted to 2D arrays, so check theinput  FITS file has 2 pixel
    #  axes (given by Nin) and 2 WCS axes (given by Nout).
    elif wcsinfo_in.Nin != 2 or wcsinfo_in.Nout != 2:
        raise ValueError("{0} is not 2-dimensional".format(hdu_in))
     
    #  Check the reference FITS file in the same way.
    elif wcsinfo_ref is None:
        raise ValueError("Failed to read WCS information from {0}".format(hdu_ref))
     
    elif wcsinfo_ref.Nin != 2 or wcsinfo_ref.Nout != 2:
        raise ValueError("{0} is not 2-dimensional".format(hdu_ref))
     
    #  Proceed if the WCS information was read OK.
 
    #  Attempt to get a mapping from pixel coords in the input FITS file to pixel
    #  coords in the reference fits file, with alignment occuring by preference in
    #  the current WCS frame. Since the pixel coordinate frame will be the base frame
    #  in each Frameset, we first invert the FrameSets. This is because the Convert method
    #  aligns current Frames, not base frames.
    wcsinfo_in.invert()
    wcsinfo_ref.invert()
    alignment_fs = wcsinfo_in.convert(wcsinfo_ref)

    #  Invert them again to put them back to their original state (i.e.
    #  base frame = pixel coords, and current Frame = WCS coords).
    wcsinfo_in.invert()
    wcsinfo_ref.invert()

    #  Check alignment was possible.
    if alignment_fs is None:
        raise Exception("Cannot find a common coordinate system shared by {0} and {1}".format(hdu_in,hdu_ref))
      
    else:
        #  Get the lower and upper bounds of the input image in pixel indices.
        #  All FITS arrays by definition have lower pixel bounds of [1,1] (unlike
        #  NDFs). Note, unlike fits AST uses FITS ordering for storing pixel axis
        #  values in an array (i.e. NAXIS1 first, NAXIS2 second, etc).
        lbnd_in = [1, 1]
        ubnd_in = [fitschan_in["NAXIS1"], fitschan_in["NAXIS2"]]
 
        #  Find the pixel bounds of the input image within the pixel coordinate
        #  system of the reference fits file.
        (lb1, ub1, xl, xu) = alignment_fs.mapbox(lbnd_in, ubnd_in, 1)
        (lb2, ub2, xl, xu) = alignment_fs.mapbox(lbnd_in, ubnd_in, 2)
 
        #  Calculate the bounds of the output image.
        lbnd_out = [int(lb1), int(lb2)]
        ubnd_out = [int(ub1), int(ub2)]
      
        #  Unlike NDFs, FITS images cannot have an arbitrary pixel origin so
        #  we need to ensure that the bottom left corner of the input image
        #  gets mapped to pixel [1,1] in the output. To do this we, extract the
        #  mapping from the alignment FrameSet and add on a ShiftMap (a mapping
        #  that just applies a shift to each axis).
        shift = [1 - lbnd_out[0],
                 1 - lbnd_out[1]]
      
        alignment_mapping = alignment_fs.getmapping()
        shiftmap = Ast.ShiftMap(shift)
        total_map = Ast.CmpMap(alignment_mapping, shiftmap)
      
        #  Modify the pixel bounds of the output image to take account of this
        #  shift of origin.
        lbnd_out[0] += shift[0]
        lbnd_out[1] += shift[1]
        ubnd_out[0] += shift[0]
        ubnd_out[1] += shift[1]
 
        #  Get the value used to represent missing pixel values
        if "BLANK" in fitschan_in:
            badval = fitschan_in["BLANK"]
            flags = Ast.USEBAD
        else:
            badval = 0
            flags = 0
 
        # Resample the data array using the above mapping.
        # total_map was pixmap; is this right?
        (npix, out, out_var) = total_map.resample(lbnd_in, ubnd_in,
                                                  hdu_in.data, None,
                                                  Ast.LINEAR, None, flags,
                                                  0.05, 1000, badval, lbnd_out,
                                                  ubnd_out, lbnd_out, ubnd_out)
 
        #  Store the aligned data in the primary HDU, and update the NAXISi keywords
        #  to hold the number of pixels along each edge of the rotated image.
        hdu_in.data = out
        fitschan_in["NAXIS1"] = ubnd_out[0] - lbnd_out[0] + 1
        fitschan_in["NAXIS2"] = ubnd_out[1] - lbnd_out[1] + 1

    #  The WCS to store in the output is the same as the reference WCS
    #  except for the extra shift of origin. So use the above shiftmap to
    #  remap the pixel coordinate frame in the reference WCS FrameSet. We
    #  can then use this FrameSet as the output FrameSet.
    wcsinfo_ref.remapframe(Ast.BASE, shiftmap)
 
    #  Attempt to write the modified WCS information to the primary HDU (i.e.
    #  convert the FrameSet to a set of FITS header cards stored in the
    #  FITS file). Indicate that we want to use original flavour of FITS-WCS.
    fitschan_in.Encoding = encoding
    fitschan_in.clear('Card')

    if fitschan_in.write(wcsinfo_ref) == 0:
        raise Exception("Failed to convert the aligned WCS to Fits-WCS")
 
    #  If successful, force the FitsChan to copy its contents into the
    #  fits header, then write the changed data and header to the output
    #  FITS file.
    else:
        fitschan_in.writefits()

    if outname is not None:
        hdu_in.writeto(outname, clobber=clobber)
    
    return hdu_in
