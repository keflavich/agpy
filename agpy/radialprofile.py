import numpy as np

def azimuthalAverage(image, center=None, stddev=False):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    Taken from http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
    via Jessica R. Lu
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

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
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr
    
    if stddev:
        # find azimuthal standard deviation
        r_int[r_int==r_int.max()] = r_int.max() - 1  # set last bin equal to second-to-last bin
        rad_mean = radial_prof[r_int]
        rad_diffsum = np.cumsum( (i_sorted-rad_mean)**2 )
        rad_std = (rad_diffsum[rind[1:]] - rad_diffsum[rind[:-1]]) / nr
        
        return rad_std

    else:
        return radial_prof
