import numpy as np

def azimuthalAverage(image, center=None, stddev=False, returnradii=False, return_nr=False ):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    stddev - if specified, return the azimuthal standard deviation instead of the average
    returnradii - if specified, return (radii_array,radial_profile)
    return_nr   - if specified, return number of pixels per radius *and* radius
    
    Taken from http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
    via Jessica R. Lu
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    deltar[-1] = 1                   # include outermost points
    rind = np.where(deltar)[0] + 1   # location of changed radius (minimum 1 point)
    nr = np.concatenate([rind[0:1], rind[1:] - rind[:-1]]) # number of radius bin
                                     # concatenate to include center pixel / innermost bin

    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = np.concatenate([csim[0:1], csim[rind[1:]] - csim[rind[:-1]]]) # include innermost bin

    radial_prof = tbin / nr

    if stddev:
        # find azimuthal standard deviation
        r_int[r_int==r_int.max()] = r_int.max() - 1  # set last bin equal to second-to-last bin
        rad_mean = radial_prof[r_int]
        rad_diffsum = np.cumsum( (i_sorted-rad_mean)**2 )
        rad_std = np.sqrt( ( np.concatenate([rad_diffsum[0:1], rad_diffsum[rind[1:]] - rad_diffsum[rind[:-1]]]) ) / nr )
        
        if returnradii: 
            return r_int[rind],rad_std
        elif return_nr:
            return nr,r_int[rind],rad_std
        else:
            return rad_std

    else:
        if returnradii: 
            return r_int[rind],radial_prof
        elif return_nr:
            return nr,r_int[rind],radial_prof
        else:
            return radial_prof
