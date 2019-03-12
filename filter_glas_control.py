#! /usr/bin/env python

#Filter preprocessed ICESat-1 GLAS points for a given input raster
#PMedit: This script has been heavily edited to accomodate our forest structure DEM Workflow
# Here, we read in:
#    a previously prepared glas csv (*asp.csv; from pc_align_prep.sh)
#    out-DEM_4m.tif that has a previously completed out-DEM_4m_control.tif (from dem_control.py)

# Slope filtering set to 'false' because it was done in dem_control.py
# Also expects out-DEM_4m_control.tif output from dem_control.py

#Then run filter_glas_forest.py for a single DEM fn or a list of DEM fn
#parallel --jobs 16 --delay 1 '~/src/demcoreg/demcoreg/filter_glas_control.py {}' ::: */dem*/out-DEM_4m.tif

import sys
import os

import numpy as np
from osgeo import gdal
from pygeotools.lib import geolib, iolib, malib, timelib

#import matplotlib.pyplot as plt
import dem_control

from imview.lib import gmtColormap, pltlib
cpt_rainbow = gmtColormap.get_rainbow()

#PMedit: site = 'hma'

#Minimum number of points required to write out _ref.csv
min_pts = 2

#Maximum value of surface slope to use
max_slope = 20.

#Reference DEM for masking
refdem_filt_list = ['/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/tandemx/TDM90/mos/TDM1_90m_circ_DEM.vrt', -15, 15]

pt_srs = geolib.wgs_srs
# Python indexing starts with 0
# This is time column: but we need to convert to YYYYMMDD
tcol = 3    # this is needed for specifying the min of the col range
xcol = 5    #lon
ycol = 4    #lat
zcol = 8    #elev_ground; also this is the max of the col range

#GLAS GLA14 Record Metadata:
#https://nsidc.org/data/docs/daac/glas_altimetry/gla14_records.html
#https://nsidc.org/data/glas/data-dictionary-glah14
#
#GLAS shot quality indicators
#----------------------------
FRir_val = 15	        #valid equal to this value:     indicates 'cloudless' waveforms
SatNdx_thresh = 2	    #valid less than this value:	indicates signals that are not 'saturated'
cld1_mswf_thresh = 15	#valid less than this value:     indicates signals not affected by multiple scattering	
wflen_thresh = 50		#valid less than this value:    indicates the total length (m) of the waveform

satndxcol = 14 
cldcol = 17
FRircol = 18
wflencol = 25

mincol = tcol
maxcol = wflencol

tcol = tcol - mincol
xcol = xcol - mincol
ycol = ycol - mincol
zcol = zcol - mincol

satndxcol = satndxcol - mincol
cldcol = cldcol - mincol
FRircol = FRircol - mincol
wflencol = FRircol - mincol


#Padding in pixels for sample radius
#PMedit: we use 4m DEM, sp pad with 8 pixels //Since we're likely dealing with 32-m products here, can just use pad=1
pad = 8
#pad = 'glas'

#PMedit:glas_dir = '/nobackupp8/deshean/icesat_glas'
glas_fn = sys.argv[1]
print(glas_fn)
glas_dir, ext = os.path.split(glas_fn)   #'/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/glas/misc/tiles_5deg_old/csv_files'
ext = os.path.splitext(ext)[0]           #'gla14_N60-70_asp'
#PMedit:ext = 'GLAH14_%s_refdemfilt' % site        

glas_npz_fn = os.path.join(glas_dir, ext+'.npz')

if not os.path.exists(glas_npz_fn):
    glas_csv_fn = os.path.splitext(glas_npz_fn)[0]+'.csv'
    print("Loading csv: %s" % glas_csv_fn)
    glas_pts = np.loadtxt(glas_csv_fn, delimiter=',', skiprows=1, dtype=None, usecols=range(mincol,maxcol+1))

    print("Saving npz: %s" % glas_npz_fn)
    np.savez_compressed(glas_npz_fn, glas_pts)

else:
    #This takes ~5 seconds to load ~9M records with 8 fields
    print("Loading npz: %s" % glas_npz_fn)
    glas_pts = np.load(glas_npz_fn)['arr_0']

dem_fn_list = sys.argv[2:]

for n,dem_fn in enumerate(dem_fn_list):

    print("%i of %i" % (n+1, len(dem_fn_list)))
    #Lat/lon extent filter
    print("Loading DEM: %s" % dem_fn)
    dem_ds = gdal.Open(dem_fn)
    dem_ma = iolib.ds_getma(dem_ds)
    dem_extent_wgs84 = geolib.ds_extent(dem_ds, t_srs=pt_srs)
    xmin, ymin, xmax, ymax = dem_extent_wgs84
    
    print("Applying spatial & quality filter") 
    x      = glas_pts[:,xcol]
    y      = glas_pts[:,ycol]
    satndx = glas_pts[:,satndxcol]
    cld    = glas_pts[:,cldcol]
    FRir   = glas_pts[:,FRircol]
    wflen  = glas_pts[:,wflencol]

    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & (satndx < SatNdx_thresh) & (cld < cld1_mswf_thresh) & (FRir == FRir_val) & (wflen < wflen_thresh))
    #idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))

    if idx.nonzero()[0].size == 0:
        print("No points after spatial & quality filtering")
        # this 'continue' makes script go to next item in for loop instead of proceeding through the rest of code below
        continue

    print("Sampling DEM at masked point locations") 
    glas_pts_fltr = glas_pts[idx]

    #print("Check rows", glas_pts_fltr[0:]) 
    print("Check row length", glas_pts_fltr.shape[1]) 

    print("Writing out %i points after spatial filter" % glas_pts_fltr.shape[0]) 
    out_csv_fn = os.path.splitext(dem_fn)[0]+'_%s.csv' % ext

    # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84 
    fmt = '%0.8f, %i, %0.6f, %0.6f, %0.2f'

    # Old gla14 csv file; before pc_align_prep.sh:
    # rec_ndx,shotn,date,lat,lon,elev,elev_geoid,elev_ground,rh100,rh50,wflen
    fmt = '%i, %i, %0.2f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f'

    if glas_pts_fltr.shape[1] == 23:    
        fmt = '%0.2f, %0.6f, %0.6f, %0.2f, %0.2f, %0.2f'
        fmt += ', %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f'
    if glas_pts_fltr.shape[1] == 7:
        # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84, z_refdem_med_WGS84, z_refdem_nmad
        fmt += ', %0.2f, %0.2f'
    elif glas_pts_fltr.shape[1] == 8:
        # dt_ordinal, dt_YYYYMMDD, lat, lon, z_WGS84, z_refdem_med_WGS84, z_refdem_nmad, lulc
        fmt += ', %0.2f, %0.2f, %i'
    elif glas_pts_fltr.shape[1] == 3:
        # After pc_align_prep.sh: lat, lon, elev_ground
        fmt = '%0.2f, %0.2f, %i'
    np.savetxt(out_csv_fn, glas_pts_fltr, fmt=fmt, delimiter=',')

    x_fltr = glas_pts_fltr[:,xcol]
    y_fltr = glas_pts_fltr[:,ycol]
    z_fltr = glas_pts_fltr[:,zcol]

    dem_mask_fn = os.path.splitext(dem_fn)[0]+'_control.tif'

    if os.path.exists(dem_mask_fn):
        print("Loading Masked DEM: %s" % dem_mask_fn)
        dem_mask_ds = gdal.Open(dem_mask_fn) 
        dem_mask = iolib.ds_getma(dem_mask_ds) 
    else:
        # Create mask here; INCOMPLETE; dont know how to specify flags (eg '--no-toamask') in call below
        #dem_mask_fn = dem_control.main(dem_fn, filt_param=refdem_filt_list)
        dem_mask_ds = dem_ds
        dem_mask = dem_ma

    #Convert input xy coordinates to raster coordinates
    mX_fltr, mY_fltr, mZ = geolib.cT_helper(x_fltr, y_fltr, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr, pY_fltr = geolib.mapToPixel(mX_fltr, mY_fltr, dem_mask_ds.GetGeoTransform())
    pX_fltr = np.atleast_1d(pX_fltr)
    pY_fltr = np.atleast_1d(pY_fltr)

    #Sample raster
    #This returns median and mad for ICESat footprint
    samp = geolib.sample(dem_mask_ds, mX_fltr, mY_fltr, pad=pad)
    samp_idx = ~(np.ma.getmaskarray(samp[:,0]))
    npts = samp_idx.nonzero()[0].size

    if npts < min_pts:
        print("Not enough points after sampling valid pixels, post control.tif mask (%i < %i)\n" % (npts, min_pts))
        continue
       
    if False:
        print("Applying slope filter, masking points with slope > %0.1f" % max_slope)
        slope_ds = geolib.gdaldem_mem_ds(dem_mask_ds, processing='slope', returnma=False)
        slope_samp = geolib.sample(slope_ds, mX_fltr, mY_fltr, pad=pad)
        slope_samp_idx = (slope_samp[:,0] <= max_slope).data
        samp_idx = np.logical_and(slope_samp_idx, samp_idx)

    npts = samp_idx.nonzero()[0].size
    if npts < min_pts:
        print("Not enough points after %0.1f deg slope mask (%i < %i)" % (max_slope, npts, min_pts))
        continue

    glas_pts_fltr_mask = glas_pts_fltr[samp_idx]
    glas_pts_fltr_mask_asp =  glas_pts_fltr_mask[:,[ycol,xcol,zcol]]

    if os.path.exists(dem_mask_fn):
        print("Writing out %i points after mask" % glas_pts_fltr_mask.shape[0]) 
        out_csv_fn_mask = os.path.splitext(out_csv_fn)[0]+'_ref.csv'
        #lat,lon,elev_ground for pc_align
        out_csv_fn_mask_asp = os.path.splitext(out_csv_fn)[0]+'_ref_asp.csv'
        #Could add DEM samp columns here
        np.savetxt(out_csv_fn_mask, glas_pts_fltr_mask, fmt=fmt, delimiter=',')
        np.savetxt(out_csv_fn_mask_asp, glas_pts_fltr_mask_asp, fmt='%0.6f, %0.6f, %0.2f', delimiter=',')

    x_fltr_mask = glas_pts_fltr_mask[:,xcol]
    y_fltr_mask = glas_pts_fltr_mask[:,ycol]
    z_fltr_mask = glas_pts_fltr_mask[:,zcol]
    mX_fltr_mask, mY_fltr_mask, mZ = geolib.cT_helper(x_fltr_mask, y_fltr_mask, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr_mask, pY_fltr_mask = geolib.mapToPixel(mX_fltr_mask, mY_fltr_mask, dem_mask_ds.GetGeoTransform())
    pX_fltr_mask = np.atleast_1d(pX_fltr_mask)
    pY_fltr_mask = np.atleast_1d(pY_fltr_mask)

    dz = z_fltr_mask - samp[samp_idx,0]

    if False:
        print "Creating plot of %i output points" % x_fltr.shape[0]
        fig_kw = {'figsize':(10,7.5)}
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True, sharey=True, **fig_kw)

        #Plot DEM color shaded relief
        hs_ma = geolib.gdaldem_wrapper(dem_fn)
        hs_clim = malib.calcperc(hs_ma, perc=(0.5, 99.5))
        dem_clim = malib.calcperc(dem_ma)
        ax1.imshow(hs_ma, cmap='gray', clim=hs_clim)
        im1 = ax1.imshow(dem_ma, cmap=cpt_rainbow, clim=dem_clim, alpha=0.5)
        cbar = pltlib.add_cbar(ax1, im1, label='DEM Elev. (m WGS84)')
       
        #Plot all color points over shaded relief
        im2 = ax2.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        #Plot all points in black
        sc2 = ax2.scatter(pX_fltr, pY_fltr, s=0.5, c='k', edgecolors='none')
        #Plot valid in color
        c = z_fltr_mask 
        sc2 = ax2.scatter(pX_fltr_mask, pY_fltr_mask, s=0.5, c=c, cmap=cpt_rainbow, vmin=dem_clim[0], vmax=dem_clim[1], edgecolors='none')
        cbar = pltlib.add_cbar(ax2, sc2, label='Pt Elev. (m WGS84)')

        #Plot time
        c = glas_pts_fltr[:,tcol]
        c_decyear = timelib.np_dt2decyear(timelib.np_o2dt(c))
        c = c_decyear
        #vmin = c.min()
        #vmax = c.max()
        vmin = 2003.14085699
        vmax = 2009.77587047
        #vmin = 20030220
        #vmax = 20091011
        im3 = ax3.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc3 = ax3.scatter(pX_fltr, pY_fltr, s=1, c=c, vmin=vmin, vmax=vmax, edgecolors='none')
        #cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year', cbar_kwargs={'format':'%0.2f'})
        cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year')

        #Plot dz
        c = dz
        vmin, vmax = malib.calcperc(c, perc=(5, 95))
        absmax = np.max(np.abs([vmin, vmax]))
        vmin = -absmax
        vmax = absmax
        im4 = ax4.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc4 = ax4.scatter(pX_fltr_mask, pY_fltr_mask, s=2, c=c, cmap='RdYlBu', vmin=vmin, vmax=vmax, edgecolors='none')
        cbar = pltlib.add_cbar(ax4, sc4, label='GCP - DEM (m)')

        for ax in (ax1, ax2, ax3, ax4):
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.set_aspect('equal', 'box-forced')

        title='%s \n %i valid points (%i initial)' % (os.path.splitext(os.path.split(dem_fn)[1])[0], pX_fltr_mask.shape[0], pX_fltr.shape[0])
        fig.suptitle(title)
        fig.tight_layout()
        #This adjusts subplots to fit suptitle
        plt.subplots_adjust(top=0.92)
        fig_fn = os.path.splitext(out_csv_fn)[0]+'.png'
        print "Saving figure: %s" % fig_fn
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close(fig)