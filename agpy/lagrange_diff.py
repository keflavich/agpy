"""
Implementation of the Lagrange differentiation scheme employed by Molinari et
al 2010 and Molinari et al 2011 to identify filamentary structures.

How do you compute the lagrange derivative coefficients analytically?  This procedure
is orders of magnitude too slow
"""

import numpy as np
from scipy.interpolate.interpolate import lagrange
import scipy.signal as signal
from agpy import timer

@timer.print_timing
def lagrange_interp(image,direction='vertical',npix=5, edge_value=0.0):

    lagrangian = np.array([[0,0,1,0,0],[0,0,-8,0,0],[0,0,0,0,0],[0,0,8,0,0],[0,0,-1,0,0]])
    if direction == 'vertical':
        pass
    elif direction == 'horizontal':
        lagrangian = lagrangian.swapaxes(0,1)
    else:
        lagrangian[:]=0
        lagrangian[0,0] = 1.0
        lagrangian[1,1] = -8.0
        lagrangian[3,3] = 8.0
        lagrangian[4,4] = -1.0


    firstderiv = signal.convolve2d(image,lagrangian,mode='same',boundary='fill')
    secondderiv = signal.convolve2d(firstderiv,lagrangian,mode='same',boundary='fill')

    return secondderiv

    yy,xx = np.indices(image.shape)
    
    pix_indices_x = np.concatenate([[xx+ii] for ii in np.arange(npix)-npix/2])
    pix_indices_y = np.concatenate([[yy]    for ii in np.arange(npix)-npix/2])
    ignore_pix_indices = (pix_indices_x < 0) + (pix_indices_x >= image.shape[0])
    pix_indices_valid_x = pix_indices_x * (True-ignore_pix_indices)
    pix_values = image[pix_indices_y,pix_indices_valid_x]
    pix_values[ignore_pix_indices] = edge_value

    
    pix_indices_flat = np.reshape(pix_indices_x,[image.shape[0]*image.shape[1],npix])
    pix_values_flat  = np.reshape(pix_values   ,[image.shape[0]*image.shape[1],npix])
    pix_indices_flat = np.array([np.arange(npix)-npix/2 for ii in xrange(image.shape[0]*image.shape[1])])

    L_polys = map(lagrange,pix_indices_flat,pix_values_flat)
    L_der2  = [L.deriv().deriv() for L in L_polys]

    pix_derivs = np.array([L(0) for L in L_der2])

    deriv_map = np.reshape(pix_derivs,image.shape)

    return deriv_map

def laplace_interp(image, splineimage=None, direction='vertical', splinesmooth=1.0):

    if splineimage is None:
        splineimage = signal.cspline2d(image,splinesmooth)
    
    laplacian = np.zeros([3,3],dtype='float32')
    laplacian[1,1] = -2.0
    if direction in ("diagonal1","diagonalleft"):
        laplacian[0,0] = 1.0
        laplacian[2,2] = 1.0
    elif direction in ("diagonal2","diagonalright"):
        laplacian[2,0] = 1.0
        laplacian[0,2] = 1.0
    elif direction == 'vertical':
        laplacian[1,0] = 1.0
        laplacian[1,2] = 1.0
    elif direction == 'horizontal':
        laplacian[0,1] = 1.0
        laplacian[2,1] = 1.0

    filtered_image = signal.convolve2d(splineimage,laplacian,mode='same',boundary='symm')

    # edges have no curvature
    filtered_image[0,:] = 0.0
    filtered_image[-1,:] = 0.0
    filtered_image[:,0] = 0.0
    filtered_image[:,-1] = 0.0

    return filtered_image

def max_curvature(image, splinesmooth=2.0):

    splineimage = signal.cspline2d(image,splinesmooth)
    curvarr = np.array([laplace_interp(image,splineimage,direction=D) for D in ['vertical','horizontal','diagonal1','diagonal2']]) * -1.0

    return curvarr.max(axis=0)
