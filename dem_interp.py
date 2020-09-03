#!/usr/bin/env python
#
# Utility to replace the nan values with interpolated values in a DEM array
# authors: B. Osmanoglu, P. Montesano

import sys
import os
import argparse
import numpy as np
from scipy import optimize

from pygeotools.lib import iolib

def frankotchellappa(dzdx,dzdy):
    '''frankotchellappa(dzdx,dzdy):
    '''
    dS = dzdx.shape
    cols = dS[1]+2   #zero pad the gradients.
    rows = dS[0]+2
    [wx, wy] = np.meshgrid((np.r_[1:cols+1]-(np.fix(cols/2)+1))/(cols-np.mod(cols,2)), (np.r_[1:rows+1]-(np.fix(rows/2)+1))/(rows-np.mod(rows,2))) 
    wx = np.fft.ifftshift(wx)
    wy = np.fft.ifftshift(wy)

    DZDX = np.fft.fft2(dzdx, (rows, cols))
    DZDY = np.fft.fft2(dzdy, (rows, cols))

    eps = np.finfo(np.double).eps

    Z = (-1j*wx*DZDX -1j*wy*DZDY)/(wx**2 + wy**2 + eps)
    z = np.fft.ifft2(Z, (rows, cols)).real
    z = z[1:-1,1:-1]

    return z

def frankotchellappaosmanoglu(dzdx,dzdy):
    '''frankotchellappaosmanoglu(dzdx,dzdy):
    '''
    z = frankotchellappa(dzdx, dzdy)
    #zr = basic.rescale(z,[-1,1])
    g1,g2 = np.gradient(z)         #g1=gy, g2=gx
    z2 = frankotchellappa(g2,g1)
    ws = (z.max()-z.min()) / (z2.max()-z2.min())
    #ws = (z.std()) / (z2.std())
    #ws = (6.*z.std()) / (2.)

    return frankotchellappa(dzdx*ws,dzdy*ws)

def fitSurface(x,y,z, weight=None, order=1):
    """planefit, fitfunc=fitSurface(x,y,z, weight=None, order=1)
    x,y,z: vectors
    weight: vector same size as x
    Known Bugs: If data contains nan values, fit might fail quietly.
    """
    if any(np.isnan(z)):
        print("Data contains nan values. Fit might fail.");
    if weight is not None:
        w = weight;
    else:
        w = np.ones(len(x))
    if order == 1:
        p0 = [0,0.1,0.1]
        fitfunc = lambda p, x, y: p[0]+p[1]*x+p[2]*y
    elif order==2:
        p0 = [0,0.1,0.1,0.1,0.1,0.1]
        fitfunc = lambda p, x, y: p[0]+p[1]*x+p[2]*y+p[3]*x**2+p[4]*y**2+p[5]*x*y
    else:
        print "order has to be 1 or 2."
        return -1

    errfunc = lambda p, x, y, z, w: abs(w*(fitfunc(p,x,y) - z))
    planefit, success = optimize.leastsq(errfunc, p0, args=(x,y,z,w))

    return planefit, fitfunc

def getparser():
    parser = argparse.ArgumentParser(description="Utility to replace the nan values with interpolated values in a DEM array")
    parser.add_argument('dem_fn', type=str, help='DEM filename')
    parser.add_argument('-weight', type=int, default=None, help='Surface fitting weight')
    parser.add_argument('-order', type=int, default=1, help='Surface fitting polynomial order')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()
 
    dem_fn = args.dem_fn
    weight = args.weight
    poly_order = args.order

    # Get original DEM
    array_dem = iolib.fn_getma(dem_fn)
    dem_ds = iolib.fn_getds(dem_fn)

    # Create empty meshgrid
    x = range (0, dem_ds.RasterXSize, 1)
    y = range (0, dem_ds.RasterYSize, 1)
    X, Y = np.meshgrid(x, y)

    print "\nGet the gradient of the dem..."
    dzdy,dzdx = np.gradient(array_dem)

    print "\nCreate copies for modification, and set the masked areas to 0 in the gradients..."
    dzdx_mod = dzdx.copy()
    dzdy_mod = dzdy.copy()
    
    print "\nOsmanoglu Algorithm..."
    #routine to get the scale corrected surface back from gradients.
    #Note that this is with the original slopes, I want the result to be close to the original DEM as possible
    fco_dem = frankotchellappaosmanoglu(dzdx, dzdy) ; #FrankotChellappa removes long wavelength trends

    print "\nSubtract the recovered DEM from the original and estimate a surface..."
    planefit, fitfunc = fitSurface(X.ravel(), Y.ravel(), (array_dem - fco_dem).ravel())

    print "\nAdd the surface to the masked gradient derived DEM..."
    dem_interp = frankotchellappaosmanoglu(dzdx_mod, dzdy_mod) + fitfunc(planefit, X, Y, weight, poly_order)


    print("\nWriting DEM with interpolated surfaces:")
    dst_fn = os.path.splitext(dem_fn)[0]+'_interp.tif'
    print(dst_fn)
    iolib.writeGTiff(dem_interp, dst_fn, dem_ds)

    return dst_fn

    # Return a numpy masked array
    return dem_interp

if __name__ == "__main__":
    main()