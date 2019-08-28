#!/usr/bin/env python
#
# Utility to get static control surfaces from a DEM
# In forested areas, DEMs that show canopy surfaces (ie, 'canopy DEMs' from low sun elev angles) will have these surfaces masked out
# Mask values in a canopy DEM using roughness and slope thresholds; above thresholds are masked out
# Mask water using TOA ortho and a min threshold
# Mask bad cloud elevs using a reference DEM with a min dz threshold
# Use before pc_align
#
# Must use correct python: source ~/anaconda3/bin/activate py2

import sys
import os
import subprocess
import glob
import argparse
import shutil

import numpy as np

from pygeotools.lib import iolib
from pygeotools.lib import malib
from pygeotools.lib import geolib
from pygeotools.lib import filtlib
from pygeotools.lib import warplib

import matplotlib
##https://stackoverflow.com/questions/37604289/tkinter-tclerror-no-display-name-and-no-display-environment-variable
matplotlib.use('Agg')
import matplotlib.pyplot, matplotlib.mlab, math
import scipy.stats

#TOA Terrain Ruggedness masked
def get_tri_mask(dem_ds, min_tri):
    print("\nApplying TRI filter (masking smooth values < %0.4f)" % min_tri)
    #dem = iolib.ds_getma(dem_ds)
    tri = geolib.gdaldem_mem_ds(dem_ds, 'TRI', returnma=True)
    tri_mask = np.ma.masked_less(tri, min_tri)
    #This should be 1 for valid surfaces, nan for removed surfaces
    tri_mask = ~(np.ma.getmaskarray(tri_mask))
    return tri_mask

#DEM roughness mask
def get_rough_mask(dem_ds, max_rough):
    print("\nApplying DEM roughness filter (masking values > %0.4f)" % max_rough)
    #dem = iolib.ds_getma(dem_ds)
    rough = geolib.gdaldem_mem_ds(dem_ds, 'Roughness', returnma=True)
    rough_mask = np.ma.masked_greater(rough, max_rough)
    #This should be 1 for valid surfaces, nan for removed surfaces
    rough_mask = ~(np.ma.getmaskarray(rough_mask))
    return rough_mask

#DEM slope mask
def get_slope_mask(dem_ds, max_slope):
    print("\nApplying DEM slope filter (masking values > %0.1f)" % max_slope)
    #dem = iolib.ds_getma(dem_ds)
    slope = geolib.gdaldem_mem_ds(dem_ds, 'slope', returnma=True)
    slope_mask = np.ma.masked_greater(slope, max_slope)
    #This should be 1 for valid surfaces, nan for removed surfaces
    slope_mask = ~(np.ma.getmaskarray(slope_mask))
    return slope_mask

#TOA reflectance mask
def get_toa_mask(toa_ds, min_toa):
    print("\nApplying TOA filter (masking values < %0.4f)" % min_toa)
    toa = iolib.ds_getma(toa_ds)
    toa_mask = np.ma.masked_less(toa, min_toa)
    #This should be 1 for valid surfaces, nan for removed surfaces
    toa_mask = ~(np.ma.getmaskarray(toa_mask))
    return toa_mask

def get_toa_fn(dem_fn):
    toa_fn = None
    dem_dir_list = os.path.split(os.path.abspath(dem_fn))[0].split(os.sep)
    import re
    #Get index of the top level pair directory containing toa (WV02_20140514_1030010031114100_1030010030896000)
    r_idx = [i for i, item in enumerate(dem_dir_list) if re.search('(_10)*(_10)*00$', item)]
    if r_idx:
        r_idx = r_idx[0]
        #Reconstruct dir
        dem_dir = (os.sep).join(dem_dir_list[0:r_idx+1])
        #Find toa.tif in top-level dir
        toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        # Check for *r100.xml here; if not exists, then break
        if not toa_fn:
            # My own version, with an edit to recognize ortho.tif, then use the 4m version of the ortho
            cmd = ['toa_calc.sh', dem_dir]
            print(cmd)
            subprocess.call(cmd)
            toa_fn = glob.glob(os.path.join(dem_dir, '*toa.tif'))
        toa_fn = toa_fn[0]
    return toa_fn

def get_min_gaus(ras_fn, sample_step=50, ncomp=3):
    # Get ma
    masked_array = iolib.fn_getma(ras_fn)
    # Sample ma
    masked_array= sample_ma(masked_array, sample_step)

    if masked_array is None:
        mean_min = 0
        stdev = 0
        print "No shift will be done. Masked array is None. Setting mean and stdv to 0."
    else:

        # Do gaussian fitting
        means, vars, weights = fit_gaus(masked_array, ncomp)

        sample_step_str = "%03d" % (sample_step)
        histo = matplotlib.pyplot.hist(masked_array.compressed(), 300, normed=True, color='gray', alpha = 0.5)
        #Write histogram
        fig_name = ras_fn.split('/')[-1].strip('.tif') + "_" + str(ncomp) + "_" + sample_step_str + '.png'
        i = 0

        out_means = []
        out_stdevs = []
        for w, m, c in zip(weights, means, vars):
            i += 1
            matplotlib.pyplot.plot(histo[1], w*scipy.stats.norm.pdf( histo[1], m, np.sqrt(c) ), linewidth=3)
            #matplotlib.pyplot.axis([min(masked_array.compressed()),max(masked_array.compressed()),0,1])
            gauss_num = 'Gaussian peak #%s' %(i)

            print 'Gaussian peak #%s (mean, stdv):  %s, %s' %(i, round(m,3), round(np.sqrt(c),3))

            out_means.append(m)
            out_stdevs.append(np.sqrt(c))

        matplotlib.pyplot.savefig(os.path.join(os.path.dirname(ras_fn),fig_name))
        matplotlib.pyplot.clf()
        print "Saved histogram fig:"
        print os.path.join(os.path.dirname(ras_fn),fig_name)

        # Find min
        mean_min = min(out_means)
        stdev = np.sqrt(vars[out_means.index(mean_min)])

    return mean_min, stdev

def fit_gaus(masked_array, ncomp):
    """
    Return the means and std devs of the ncomp number of gaussians computed for the histogram of an input masked array
        write out a figure of the histogram and the gaussians based on the raster filename from which the masked array is derived.
    """
    import sklearn.mixture
    # http://stackoverflow.com/questions/10143905/python-two-curve-gaussian-fitting-with-non-linear-least-squares/19182915#19182915
    X_compress = masked_array.compressed()
    X_reshape = np.reshape(X_compress,(masked_array.compressed().size,1))

    clf = sklearn.mixture.GaussianMixture(n_components=ncomp, covariance_type='full')
    clf.fit(X_reshape)

    ml = clf.means_
    wl = clf.weights_
    cl = clf.covariances_
    ms = [m[0] for m in ml]
    cs = [np.sqrt(c[0][0]) for c in cl]
    ws = [w for w in wl]

    return ms, cs, ws


def sample_ma(array, sampleStep, min_val=-99):
    """
    Return a sub-sampled masked array for histogram and guassian analysis
        Do histogram of image by regularly sampling a 'pct' of the input image's pixels
        Provides an even sample from across the entire image without having to analyze the entire array
    """
    #Creating data range
    masked_array = np.ma.masked_less_equal(array, min_val)         # mask all values inside this interval
    masked_array = np.ma.masked_invalid(masked_array)              # mask all nan and inf values

    # Numpy slicing to sample image for histogram generation
    # Get size
    nrow,ncol = masked_array.shape
    #print '\nRaster histogram: sampling & estimating gaussian peaks'
    #print 'Array dims: ' + str(nrow) + " , " + str(ncol)

    # [start:stop:step]
    print 'Sampling the rows, cols with sample step: %s' %(sampleStep)

    masked_array = masked_array[0::sampleStep,0::sampleStep]
    sz = masked_array.size
    #print '\tNum. elements in NEW sampled array: %s' %(sz)
    #print "\t: min, max, med, mean, std"
    #print "\t:",masked_array.min(),masked_array.max(),np.ma.median(masked_array),masked_array.mean(),masked_array.std()

    if masked_array.compressed().size > 1:
        return masked_array

def force_symlink(file1, file2):
    import errno
    try:
        os.symlink(file1, file2)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)

def getparser():
    parser = argparse.ArgumentParser(description="Utility to get static control surfaces from a DEM")
    parser.add_argument('dem_fn', type=str, help='DEM filename')
    parser.add_argument('-out_dir', type=str, default=None, help='Output directory')
    parser.add_argument('-ndv', default=0, type=int, help='Output nodata value (default: %(default)s)')
    parser.add_argument('-max_slope', type=int, default=20, help='Max slope (degrees) that will be included (default: %(default)s)')
    parser.add_argument('-max_rough', type=float, default=0.75, help='Max roughness (diff *input units* from adjacent) to be included (default: %(default)s)')
    parser.add_argument('-min_toa', type=float, default=0.15, help='Min TOA to be included (default: %(default)s)')
    parser.add_argument('-min_toatri', type=float, default=0.001, help='Min TOA TRI (the smoothest surfaces) to be included (default: %(default)s)')
    parser.add_argument('-filt_param', nargs=3, default=None, help='Filter parameter list (e.g., size, min max, ref_fn min max) (default: %(default)s)')
    #parser.add_argument('-dem_coarscomp_fn', type=str, default=None, help='The filename of the coarse companion to an input DEM (default: %(default)s)')
    parser.add_argument('--dilate_con', type=int, default=None, help='Dilate control mask with this many iterations (default: %(default)s)')
    parser.add_argument('--dilate_toa', type=int, default=None, help='Dilate TOA mask with this many iterations (default: %(default)s)')
    parser.add_argument('--no-filtdz', dest='filtdz', action='store_false', help='Turn dz filtering off')
    parser.set_defaults(filtdz=True)
    parser.add_argument('--no-roughmask', dest='roughmask', action='store_false', help='Turn roughness masking off')
    parser.set_defaults(roughmask=True)
    parser.add_argument('--no-toamask', dest='toamask', action='store_false', help='Turn TOA dark masking off')
    parser.set_defaults(toamask=True)
    parser.add_argument('--no-toatrimask', dest='toatrimask', action='store_false', help='Turn TOA TRI (smoothness) masking off')
    parser.set_defaults(toatrimask=True)
    parser.add_argument('--no-slopemask', dest='slopemask', action='store_false', help='Turn slope masking off')
    parser.set_defaults(slopemask=True)
    #parser.add_argument('--no-slopemask-coarse', dest='slopemaskcoarse', action='store_false', help='Turn coarse slope masking off')
    #parser.set_defaults(slopemaskcoarse=True)
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    dem_fn = args.dem_fn

    # Write out because they will be used to mask CHM
    writeall = True

    # Auto compute min TOA with gaussian mixture model
    compute_min_toa = True

    dirname, demname = os.path.split(dem_fn)

    # The subdir in which the DEM.tif sits will be the pairname
    pairname = os.path.split(dirname)[1]
    print("Pairname:", pairname)

    if args.out_dir is not None:

        # Create symlink in out_dir to: (1) original out-DEM_4m (2) *_ortho_4m.tif (3) All *.xml files
        # This should look like <out_dir>/<pairname>_out-DEM_4m
        dem_fn_lnk = os.path.join(args.out_dir, pairname + '_' +  demname)
        force_symlink(dem_fn, dem_fn_lnk)
        force_symlink(os.path.join(dirname, pairname + '_ortho_4m.tif'), os.path.join(args.out_dir, pairname + '_ortho_4m.tif') )
        xml_list = [f for f in os.listdir(dirname) if f.endswith('r100.xml')]

        print("\nSymlinks made for:")
        for x in xml_list:
            print(x)
            shutil.copy2(os.path.join(dirname,x), args.out_dir)

        out_fn_base = os.path.splitext(dem_fn_lnk)[0]

        dem_fn = dem_fn_lnk
    else:
        out_fn_base = os.path.splitext(dem_fn)[0]

    print("\nBasename for output files:")
    print(out_fn_base)

    #Need some checks on these
    param = args.filt_param
    if param is not None and len(param) == 1:
        param = param[0]
    # Get original DEM
    dem = iolib.fn_getma(dem_fn)

    print("\nLoading input DEM into masked array")
    dem_ds = iolib.fn_getds(dem_fn)

    toa_mask = None
    toa_tri_mask = None # probably not used by itself; done as part of toa_mask
    rough_mask = None
    slope_mask = None
    mask_list = [toa_tri_mask, toa_mask, rough_mask, slope_mask]

    if args.filtdz:
        print("\nFilter with dz from ref DEM to remove cloud returns and blunders (shadows)...")
        print("Reference DEM: %s" % os.path.split(param[0])[1] )
        print("Absolute dz (+/-): %s \n" % param[2] )
        #May need to cast input ma as float32 so np.nan filling works
        dem = dem.astype(np.float32)

        #Difference filter, need to specify ref_fn and range
        #Could let the user compute their own dz, then just run a standard range or absrange filter
        ref_fn = param[0]
        ref_ds = warplib.memwarp_multi_fn([ref_fn,], res=dem_ds, extent=dem_ds, t_srs=dem_ds)[0]
        ref = iolib.ds_getma(ref_ds)
        param = map(float, param[1:])

        # A dem that has been masked based on the dz filter
        dem = filtlib.dz_fltr_ma(dem, ref, rangelim=param)

        if writeall:
            #out_fn = os.path.splitext(dem_fn)[0]+'_dzfilt.tif'
            out_fn = os.path.join(out_fn_base +'_dzfilt.tif')
            print("Writing out %s\n" % out_fn)
            iolib.writeGTiff(dem, out_fn, src_ds=dem_ds, ndv=args.ndv)

    #Initialize a control mask that we'll update
    #True (1) represents "valid" unmasked pixel, False (0) represents "invalid" pixel to be masked
    controlmask = ~(np.ma.getmaskarray(dem))

    # DEM masking: Each block returns a masked output (not a mask)
    #    TOA: mask dark and/or smooth areas (shadows and/or water)
    #    Roughness
    #    Slope

    if args.toamask or args.toatrimask:
        #try:
            print("\nCompute TOA from ortho...\n")
            toa_fn = get_toa_fn(out_fn_base + '.tif') ##--->dem_fn
            print("\nWarp TOA to DEM...\n")
            print(toa_fn)
            toa_ds = warplib.memwarp_multi_fn([toa_fn,], res=dem_ds, extent=dem_ds, t_srs=dem_ds)[0]

            if args.toamask:

                if compute_min_toa:

                    # Compute a good min TOA value
                    m,s = get_min_gaus(toa_fn, 50, 4)
                    min_toa = m + s
                    min_toa = m
                else:
                    min_toa = args.min_toa

                with open(os.path.join(os.path.split(toa_fn)[0], "min_toa_" + pairname + ".txt"), "w") as text_file:
                    text_file.write(os.path.basename(__file__))
                    text_file.write("\nMinimum TOA used for mask:\n{0}".format(min_toa))

                # Should mask dark areas and dilate
                toa_mask = get_toa_mask(toa_ds, min_toa)

                #Dilate the mask
                if args.dilate_toa is not None:
                    niter = args.dilate_toa
                    print("Dilating TOA mask with %i iterations" % niter)
                    from scipy import ndimage
                    toa_mask = ~(ndimage.morphology.binary_dilation(~toa_mask, iterations=niter))

                controlmask = np.logical_and(toa_mask, controlmask)

                # Mask islands here
                controlmask = malib.mask_islands(controlmask, 5)

                if writeall:
                    #out_fn = out_fn_base+'_toamask.tif'
                    out_fn = os.path.join(out_fn_base +'_toamask.tif')
                    print("Writing out %s\n" % out_fn)
                    iolib.writeGTiff(toa_mask, out_fn, src_ds=dem_ds)

            if args.toatrimask:
                # Should mask smooth areas (measures local variance)
                toa_tri_mask = get_tri_mask(toa_ds, args.min_toatri)
                controlmask = np.logical_and(toa_tri_mask, controlmask)

                if writeall:
                    #out_fn = out_fn_base+'_toatrimask.tif'
                    out_fn = os.path.join(out_fn_base +'_toatrimask.tif')
                    print("Writing out %s\n" % out_fn)
                    iolib.writeGTiff(toa_tri_mask, out_fn, src_ds=dem_ds)

        #except Exception, e:
            #print "\tFailed to apply TOA masking.\n"

    if args.slopemask:
        slope_mask = get_slope_mask(dem_ds, args.max_slope)
        controlmask = np.logical_and(slope_mask, controlmask)

        #if args.slopemaskcoarse:

            #dem_fn2 = args.dem_coarscomp_fn

            #print("\nLoading input coarse DEM into masked array")
            #dem2_ds = iolib.fn_getds(dem_fn2)
            #slope_mask = get_slope_mask(dem2_ds, args.max_slope)
            #controlmask = np.logical_and(slope_mask, controlmask)

        if writeall:
            #out_fn = out_fn_base+'_slopemask.tif'
            out_fn = os.path.join(out_fn_base +'_slopemask.tif')
            print("Writing out %s\n" % out_fn)
            iolib.writeGTiff(slope_mask, out_fn, src_ds=dem_ds)

    # CHM mask will be a subset of the Control mask; slope_mask, toa_mask, toa_tri_mask
    chmmask = controlmask
    print("Generating final CHM mask to apply later")
    #out_fn = out_fn_base+'_chmmask.tif'
    out_fn = os.path.join(out_fn_base +'_chmmask.tif')
    print("Writing out %s\n" % out_fn)
    iolib.writeGTiff(chmmask, out_fn, src_ds=dem_ds)

    if args.roughmask:
        rough_mask = get_rough_mask(dem_ds, args.max_rough)
        controlmask = np.logical_and(rough_mask, controlmask)
        if writeall:
            #out_fn = out_fn_base+'_controlmask.tif'
            out_fn = os.path.join(out_fn_base +'_controlmask.tif')
            print("Writing out %s\n" % out_fn)
            iolib.writeGTiff(rough_mask, out_fn, src_ds=dem_ds)

    print("Generating final mask to use for reference surfaces, and applying to input DEM")
    #Now invert to use to create final masked array
    controlmask = ~controlmask

    #Dilate the mask
    if args.dilate_con is not None:
        niter = args.dilate_con
        print("Dilating control mask with %i iterations" % niter)
        from scipy import ndimage
        controlmask = ~(ndimage.morphology.binary_dilation(~controlmask, iterations=niter))

    #Apply mask to original DEM - use these surfaces for co-registration
    newdem = np.ma.array(dem, mask=controlmask)

    if True:
        print("\nStats of valid DEM with maskes applied:")
        valid_stats = malib.print_stats(newdem)
        valid_stats_med = valid_stats[5]

    print("\nWriting DEM control surfaces:")
    #if args.out_dir is not None:
    #    dst_fn = os.path.join(args.out_dir, os.path.split(dirname)[1] + os.path.splitext(demname)[0]+'_control.tif')
    #else:
    #    dst_fn = os.path.splitext(dem_fn)[0]+'_control.tif'
    dst_fn = os.path.join(out_fn_base +'_control.tif')
    print(dst_fn)
    iolib.writeGTiff(newdem, dst_fn, dem_ds)

    return dst_fn

if __name__ == "__main__":
    main()