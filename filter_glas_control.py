#! /usr/bin/env python

#Filter preprocessed ICESat-1 GLAS points for a given input raster
#
# required: source ~/anaconda3/bin/activate py2
#
# We read in:
#    a previously prepared glas csv
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

import matplotlib.pyplot as plt
import dem_control

from imview.lib import gmtColormap, pltlib
cpt_rainbow = gmtColormap.get_rainbow()

def glas_qfilt(glas_pts, satndxcol=14, cldcol=17, FRircol=18, wflencol=25):
    """Do quality filtering of ICESat-GLAS using standard filtering
    """
    #GLAS GLA14 Record Metadata:
    #https://nsidc.org/data/docs/daac/glas_altimetry/gla14_records.html
    #https://nsidc.org/data/glas/data-dictionary-glah14
    #
    #GLAS shot quality indicators
    #----------------------------
    FRir_val = 15	        #valid equal to this value:     indicates 'cloudless' waveforms
    SatNdx_thresh = 2	    #valid less than this value:	indicates signals that are not 'saturated'
    cld1_mswf_thresh = 15	#valid less than this value:    indicates signals not affected by multiple scattering
    wflen_thresh = 50		#valid less than this value:    indicates the total length (m) of the waveform

    print("Applying quality filter")
    satndx = glas_pts[:,satndxcol]
    cld    = glas_pts[:,cldcol]
    FRir   = glas_pts[:,FRircol]
    wflen  = glas_pts[:,wflencol]

    idx = ((satndx < SatNdx_thresh) & (cld < cld1_mswf_thresh) & (FRir == FRir_val) ) # & (wflen < wflen_thresh))
    glas_pts = glas_pts[idx]

    return glas_pts

def get_raster_idx(x_vect, y_vect, pt_srs, ras_ds, max_slope=20):
    """Get raster index corresponding to the set of X,Y locations
    """
    #Convert input xy coordinates to raster coordinates
    mX_fltr, mY_fltr, mZ = geolib.cT_helper(x_vect, y_vect, 0, pt_srs, geolib.get_ds_srs(ras_ds))
    pX_fltr, pY_fltr = geolib.mapToPixel(mX_fltr, mY_fltr, ras_ds.GetGeoTransform())
    pX_fltr = np.atleast_1d(pX_fltr)
    pY_fltr = np.atleast_1d(pY_fltr)

    #Sample raster
    #This returns median and mad for ICESat footprint
    samp = geolib.sample(ras_ds, mX_fltr, mY_fltr, pad=pad)
    samp_idx = ~(np.ma.getmaskarray(samp[:,0]))
    npts = samp_idx.nonzero()[0].size

    if False:
        print("Applying slope filter, masking points with slope > %0.1f" % max_slope)
        slope_ds = geolib.gdaldem_mem_ds(ras_ds, processing='slope', returnma=False)
        slope_samp = geolib.sample(slope_ds, mX_fltr, mY_fltr, pad=pad)
        slope_samp_idx = (slope_samp[:,0] <= max_slope).data
        samp_idx = np.logical_and(slope_samp_idx, samp_idx)

    return samp, samp_idx, npts, pX_fltr, pY_fltr

#Minimum number of points required to write out _ref.csv
min_pts = 2

#Maximum value of surface slope to use
max_slope = 20

#Reference DEM for masking
refdem_filt_list = ['/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/tandemx/TDM90/mos/TDM1_90m_circ_DEM.vrt', -15, 15]

pt_srs = geolib.wgs_srs

# Header & format from $NOBACKUP/userfs02/data/glas/circ_boreal
hdr_pre="rec_ndx,shotn,datatake,"
fmt_pre='%i, %i, %s,'

hdr="date,lat,lon,elev,elev_geoid,elev_ground,"
fmt='%0.2f, %0.6f, %0.6f, %0.3f, %0.2f, %0.2f,'
hdr+="g14_rh100,g14rh50,g14_wflen,BeamCoelv,DEM_hires_src,satNdx,satElevCorr,"
fmt+='%0.4f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f,'
hdr+="ElvuseFlg,cld1_mswf,FRir_qaFlag,lead,trail,MedH,MeanH,QMCH,Centroid,wflen,"
fmt+='%0.5f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f,'

hdr+="ht1,ht2,ht3,ht4,fslope,eratio,Senergy,"
fmt+='%0.6f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f,'
hdr+="h10,h20,h25,h30,h40,h50,h60,h70,h75,h80,h90,h100"
fmt+='%0.7f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.8f, %0.2f'

hdr_full=hdr_pre + hdr
fmt_full=fmt_pre + fmt

hdr_full_list = hdr_full.split(',')
fmt_full_list = fmt_full.split(',')

# For choosing the range of cols for loadtxt: the min max of the intial input
mincol = hdr_full_list.index('date')        # this is needed for specifying the min of the col range
maxcol = hdr_full_list.index('wflen')       # this is needed for specifying the max of the col range

# Filter out some problematic cols: some cols in the range will have rare problematic rows (eg; 2 cols combined)
hdr_filt_list = hdr_full_list[mincol:maxcol+1]

# Filter out some problematic cols: some cols in the range will have rare problematic rows (eg; 2 cols combined)
#exclude_list = [hdr_full_list.index('rec_ndx'), hdr_full_list.index('shotn'), hdr_full_list.index('datatake'), hdr_full_list.index('eratio') , hdr_full_list.index('Senergy') ]
#hdr_idx_list = [item for item in range(0,len(hdr_full_list)) if item not in exclude_list]
#hdr_idxexcl_list = [item for item in range(0,len(hdr_full_list)) if item in exclude_list]
#hdr_excl_list = [hdr_full_list[i] for i in hdr_idxexcl_list]
#hdr_filt_list = [hdr_full_list[i] for i in hdr_idx_list]
#fmt_filt_list = [fmt_full_list[i] for i in hdr_idx_list]
#print("Hdr filt list 1: ", hdr_filt_list)

# Cut down hdr and fmt lists accordingly
hdr_filt_list = hdr_full_list[mincol:maxcol + 1]
fmt_filt_list = fmt_full_list[mincol:maxcol + 1]
#print("Hdr filt list: ", hdr_filt_list)

tcol = hdr_filt_list.index('date')
xcol = hdr_filt_list.index('lon')         #lon
ycol = hdr_filt_list.index('lat')         #lat
zcol = hdr_filt_list.index('elev_ground') #elev_ground; also this is the max of the col range
satndxcol = hdr_filt_list.index('satNdx')
cldcol = hdr_filt_list.index('cld1_mswf')
FRircol = hdr_filt_list.index('FRir_qaFlag')
wflencol = hdr_filt_list.index('wflen')

hdr_out = ','.join(hdr_filt_list)
fmt_out = ','.join(fmt_filt_list)

#Padding in pixels for sample radius
#PMedit: we use 4m DEM, sp pad with 8 pixels //Since we're likely dealing with 32-m products here, can just use pad=1
pad = 8

glas_fn = sys.argv[1]
print(glas_fn)

glas_dir, ext = os.path.split(glas_fn)   #'/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/glas/misc/tiles_5deg_old/csv_files'
ext = os.path.splitext(ext)[0]           #'gla14_N60-70_asp'


glas_npz_fn = os.path.join(glas_dir, ext + '.npz')
glas_npz_qfilt_fn = os.path.join(glas_dir, ext + '_qfilt.npz')
glas_qfilt_csv =    os.path.join(glas_dir, ext + '_qfilt.csv')
print("")

if not os.path.exists(glas_npz_fn):
    # Make npz
    glas_csv_fn = os.path.splitext(glas_npz_fn)[0]+'.csv'
    print("Loading csv: %s" % glas_csv_fn)

    print("Full hdr length is %s" % (len(hdr_full_list) ) )
    print("Full hdr is: %s" % (hdr_full_list ) )

    print("Filtered hdr length is %s" % (len(hdr_filt_list) ) )
    #print("Filtered hdr removed cols: %s" % (hdr_excl_list) )

    glas_pts = np.loadtxt(glas_csv_fn, delimiter=',', skiprows=1, dtype=None, usecols=range(mincol,maxcol+1)) #hdr_idx_list)

    print("Check # cols of incoming set of GLAS", glas_pts.shape[1])
    print("Saving npz: %s" % glas_npz_fn)
    np.savez_compressed(glas_npz_fn, glas_pts)

    print("Doing quality filtering of ICESat-GLAS using standard filtering")
    glas_pts = glas_qfilt(glas_pts, satndxcol=satndxcol, cldcol=cldcol, FRircol=FRircol, wflencol=wflencol)

    print("Saving quality filtered npz: %s" % glas_npz_qfilt_fn)
    np.savez_compressed(glas_npz_qfilt_fn, glas_pts)

else:
    #Load npz
    if not os.path.exists(glas_npz_qfilt_fn):
        #Not yet quality filtered
        #This takes ~5 seconds to load ~9M records with 8 fields
        print("Loading npz: %s" % glas_npz_fn)
        glas_pts = np.load(glas_npz_fn)['arr_0']
        print("Check # cols of incoming set of GLAS", glas_pts.shape[1])

        print("Doing quality filtering of ICESat-GLAS using standard filtering")
        glas_pts = glas_qfilt(glas_pts, satndxcol=satndxcol, cldcol=cldcol, FRircol=FRircol, wflencol=wflencol)

        print("Saving quality filtered npz: %s" % glas_npz_qfilt_fn)
        np.savez_compressed(glas_npz_qfilt_fn, glas_pts)
    else:
        #Already quality filtered
        print("Loading quality filtered npz: %s" % glas_npz_qfilt_fn)
        glas_pts = np.load(glas_npz_qfilt_fn)['arr_0']


if not os.path.exists(glas_qfilt_csv):
    print("Saving quality filtered csv: %s" % glas_qfilt_csv)
    np.savetxt(glas_qfilt_csv, glas_pts, header=hdr_out, fmt=fmt_out, delimiter=',', comments='')

dem_fn_list = sys.argv[2:]

print("Check # rows of incoming qfilt set of GLAS", glas_pts.shape[0])
print("Check # cols of incoming qfilt set of GLAS", glas_pts.shape[1])

for n,dem_fn in enumerate(dem_fn_list):

    main_dir, pairname = os.path.split(os.path.split(dem_fn)[0])

    print("%i of %i" % (n+1, len(dem_fn_list)))
    #Lat/lon extent filter
    print("Loading DEM: %s" % dem_fn)
    dem_ds = gdal.Open(dem_fn)
    dem_ma = iolib.ds_getma(dem_ds)
    dem_extent_wgs84 = geolib.ds_extent(dem_ds, t_srs=pt_srs)
    xmin, ymin, xmax, ymax = dem_extent_wgs84

    print("Applying spatial filter")
    x      = glas_pts[:,xcol]
    y      = glas_pts[:,ycol]
    #satndx = glas_pts[:,satndxcol]
    #cld    = glas_pts[:,cldcol]
    #FRir   = glas_pts[:,FRircol]
    #wflen  = glas_pts[:,wflencol]

    #idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & (satndx < SatNdx_thresh) & (cld < cld1_mswf_thresh) & (FRir == FRir_val) & (wflen < wflen_thresh))
    idx = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))
    print("DEM extent (xmin,ymin,xmax,ymax): ", dem_extent_wgs84)
    print("GLAS extent (xmin,ymin,xmax,ymax): ", min(x),min(y),max(x),max(y) )
    if idx.nonzero()[0].size == 0:
        print("No points after spatial & quality filtering")
        # this 'continue' makes script go to next item in for loop instead of proceeding through the rest of code below
        continue

    print("Sampling DEM at masked point locations")
    glas_pts_fltr = glas_pts[idx]

    print("Check # rows", glas_pts_fltr.shape[0])
    print("Check # cols", glas_pts_fltr.shape[1])

    print("Writing out %i of ICESat-GLAS shots after spatial filter" % glas_pts_fltr.shape[0])
    out_csv_fn = os.path.splitext(dem_fn)[0]+'_%s.csv' % ext

    print("Writing out CSV of spatial filtered with # cols = %i" % glas_pts_fltr.shape[1] )
    np.savetxt(out_csv_fn, glas_pts_fltr, header=hdr_out, fmt=fmt_out, delimiter=',', comments='')

    x_fltr = glas_pts_fltr[:,xcol]
    y_fltr = glas_pts_fltr[:,ycol]
    z_fltr = glas_pts_fltr[:,zcol]

    # The mask used to find co-reg points (it excludes forests)
    dem_mask_fn = os.path.splitext(dem_fn)[0]+'_control.tif'
    # The mask used to sample (eg, with zonal stats) the DEM with qfilt GLAS
    dem_chmmask_fn = os.path.splitext(dem_fn)[0]+'_chmmask.tif'


    if os.path.exists(dem_mask_fn):
        print("Loading 'control' DEM (masked for co-reg): %s" % dem_mask_fn)
        dem_mask_ds = gdal.Open(dem_mask_fn)
        dem_mask = iolib.ds_getma(dem_mask_ds)
        print("Loading 'chm mask' DEM (masked for zonal stats of all valid surfaces): %s" % dem_chmmask_fn)
        dem_chmmask_ds = gdal.Open(dem_chmmask_fn)
    else:
        # Create mask here; INCOMPLETE; dont know how to specify flags (eg '--no-toamask') in call below
        #dem_mask_fn = dem_control.main(dem_fn, filt_param=refdem_filt_list)
        dem_mask_ds = dem_ds
        dem_mask = dem_ma

    # Run function to get index of raster pixels that match quality-filtered ICESat-GLAS
    coreg_samp, coreg_samp_idx, npts, pX_fltr, pY_fltr = get_raster_idx(x_fltr, y_fltr, pt_srs, dem_mask_ds, max_slope=max_slope)
    chm_samp, chm_samp_idx, npts_chm, pX_fltr_chm, pY_fltr_chm = get_raster_idx(x_fltr, y_fltr, pt_srs, dem_chmmask_ds, max_slope=max_slope)

    if npts < min_pts:
        print("Not enough points after sampling valid pixels, post control.tif mask (%i < %i)\n" % (npts, min_pts))
        continue

    #npts = coreg_samp_idx.nonzero()[0].size
    #if npts < min_pts:
    #    print("Not enough points after %0.1f deg slope mask (%i < %i)" % (max_slope, npts, min_pts))
    #    continue

    # Ground Surface (Co-reg) Set: provides the set of qfilt ICESat-GLAS used to co-reg (good GLAS that match valid pixels of ground surfaces)
    # All usecols
    glas_pts_fltr_coreg = glas_pts_fltr[coreg_samp_idx]
    # 3 cols for pc_align: x,y, z
    glas_pts_fltr_coreg_asp =  glas_pts_fltr_coreg[:,[ycol,xcol,zcol]]

    # Valid Surface Set: valid pixels of any ground and canopy surface
    glas_pts_fltr_valsurf = glas_pts_fltr[chm_samp_idx]
    print("Check # of ICESat-GLAS shots for valid surfaces", glas_pts_fltr_valsurf.shape[0])
    print("Writing out %i ICESat-GLAS shots for valid surfaces" % glas_pts_fltr_valsurf.shape[0])
    out_csv_valsurf_fn = os.path.splitext(dem_fn)[0]+'_%s_validsurfs.csv' % ext
    print("Writing out CSV of of ICESat-GLAS shots for valid surfaces with # cols = %i" % glas_pts_fltr_valsurf.shape[1] )
    np.savetxt(out_csv_valsurf_fn, glas_pts_fltr_valsurf, header=hdr_out, fmt=fmt_out, delimiter=',', comments='')

    if os.path.exists(dem_mask_fn):
        print("Writing out %i ICESat-GLAS shots for co-registration" % glas_pts_fltr_coreg.shape[0])
        out_csv_fn_coreg = os.path.splitext(out_csv_fn)[0]+'_ref.csv'
        #lat,lon,elev_ground for pc_align
        out_csv_fn_coreg_asp = os.path.splitext(out_csv_fn)[0]+'_ref_asp.csv'
        #Could add DEM samp columns here
        np.savetxt(out_csv_fn_coreg, glas_pts_fltr_coreg, header=hdr_out, fmt=fmt_out, delimiter=',', comments='')
        np.savetxt(out_csv_fn_coreg_asp, glas_pts_fltr_coreg_asp, fmt='%0.6f, %0.6f, %0.2f', delimiter=',')

    # For plotting the qfilt ICESat-GLAS used for co-reg
    x_fltr_mask_coreg = glas_pts_fltr_coreg[:,xcol]
    y_fltr_mask_coreg = glas_pts_fltr_coreg[:,ycol]
    z_fltr_mask_coreg = glas_pts_fltr_coreg[:,zcol]
    mX_fltr_mask_coreg, mY_fltr_mask_coreg, mZ_coreg = geolib.cT_helper(x_fltr_mask_coreg, y_fltr_mask_coreg, 0, pt_srs, geolib.get_ds_srs(dem_mask_ds))
    pX_fltr_mask_coreg, pY_fltr_mask_coreg = geolib.mapToPixel(mX_fltr_mask_coreg, mY_fltr_mask_coreg, dem_mask_ds.GetGeoTransform())
    pX_fltr_mask_coreg = np.atleast_1d(pX_fltr_mask_coreg)
    pY_fltr_mask_coreg = np.atleast_1d(pY_fltr_mask_coreg)

    # For plotting the qfilt ICESat-GLAS used for examining all valid surfaces
    x_fltr_mask_valsurf = glas_pts_fltr_valsurf[:,xcol]
    y_fltr_mask_valsurf = glas_pts_fltr_valsurf[:,ycol]
    z_fltr_mask_valsurf = glas_pts_fltr_valsurf[:,zcol]
    mX_fltr_mask_valsurf, mY_fltr_mask_valsurf, mZ_valsurf = geolib.cT_helper(x_fltr_mask_valsurf, y_fltr_mask_valsurf, 0, pt_srs, geolib.get_ds_srs(dem_chmmask_ds))
    pX_fltr_mask_valsurf, pY_fltr_mask_valsurf = geolib.mapToPixel(mX_fltr_mask_valsurf, mY_fltr_mask_valsurf, dem_chmmask_ds.GetGeoTransform())
    pX_fltr_mask_valsurf = np.atleast_1d(pX_fltr_mask_valsurf)
    pY_fltr_mask_valsurf = np.atleast_1d(pY_fltr_mask_valsurf)

    # Get the elev dif b/w GLAS and DEM
    dz = z_fltr_mask_coreg - coreg_samp[coreg_samp_idx,0]

    if True:
        print "Creating plot of %i output points" % x_fltr.shape[0]
        fig_kw = {'figsize':(10,7.5)}
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharex=True, sharey=True, **fig_kw)

        #Plot DEM color shaded relief
        #
        hs_ma = geolib.gdaldem_wrapper(dem_fn)
        hs_clim = malib.calcperc(hs_ma, perc=(0.5, 99.5))
        dem_clim = malib.calcperc(dem_ma)
        ax1.imshow(hs_ma, cmap='gray', clim=hs_clim)
        im1 = ax1.imshow(dem_ma, cmap=cpt_rainbow, clim=dem_clim, alpha=0.5)
        cbar = pltlib.add_cbar(ax1, im1, label='DEM Elev. (m WGS84)')

        #Plot all points in extent over shaded relief; overplot coreg points
        #
        ax2.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        # Show 4m ortho
        toa_fn = dem_control.get_toa_fn(dem_fn)
        print("Loading TOA ortho: %s" % toa_fn)
        toa_ds = gdal.Open(toa_fn)
        toa_ma = iolib.ds_getma(toa_ds)
        toa_clim = malib.calcperc(toa_ma, perc=(0.5, 99.5))
        im2 = ax2.imshow(toa_ma, cmap='gray', clim=toa_clim, alpha=1)
        # Show roughmask (should be called 'controlmask' which shows control surfaces
        #dem_roughmask_fn = os.path.splitext(dem_fn)[0]+'_roughmask.tif'
        #dem_roughmask_ds = gdal.Open(dem_roughmask_fn)
        #dem_roughmask_ma = iolib.ds_getma(dem_roughmask_ds)
        #dem_chmmask_ma = iolib.ds_getma(dem_chmmask_ds)
        #im2 = ax2.imshow(dem_chmmask_ma, cmap='Dark2', clim=malib.calcperc(dem_chmmask_ma), alpha=0.5)
        #im2 = ax2.imshow(dem_roughmask_ma, cmap='gray', clim=malib.calcperc(dem_roughmask_ma), alpha=1)
        #Plot all points in black
        sc2 = ax2.scatter(pX_fltr, pY_fltr, s=0.25, c='k', edgecolors='#66bd63')
        #Plot valid for co-reg in color
        ##sc2 = ax2.scatter(pX_fltr_mask_coreg, pY_fltr_mask_coreg, s=3, c=z_fltr_mask_coreg, cmap=cpt_rainbow, vmin=dem_clim[0], vmax=dem_clim[1], edgecolors='none')
        ##cbar = pltlib.add_cbar(ax2, sc2, label='ICESat-GLAS Elev. (m WGS84)')

        #Plot valid surface points over shaded relief; overplot coreg points
        #
        im3 = ax3.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        #Plot valid surface points in black
        sc3 = ax3.scatter(pX_fltr_mask_valsurf, pY_fltr_mask_valsurf, s=0.5, c='k', edgecolors='none')
        #Plot valid for co-reg in color
        sc3 = ax3.scatter(pX_fltr_mask_coreg, pY_fltr_mask_coreg, s=3, c=z_fltr_mask_coreg, cmap=cpt_rainbow, vmin=dem_clim[0], vmax=dem_clim[1], edgecolors='none')
        cbar = pltlib.add_cbar(ax3, sc3, label='ICESat-GLAS Elev. (m WGS84)')

        #Plot time
        ##c = glas_pts_fltr[:,tcol]
        ##c_decyear = timelib.np_dt2decyear(timelib.np_o2dt(c))
        ##c = c_decyear
        #vmin = c.min()
        #vmax = c.max()
        ##vmin = 2003.14085699
        ##vmax = 2009.77587047
        ##im3 = ax3.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        ##sc3 = ax3.scatter(pX_fltr_mask_valsurf, pY_fltr_mask_valsurf, s=1, c=z_fltr_mask_valsurf, vmin=vmin, vmax=vmax, edgecolors='none')
        ##cbar = pltlib.add_cbar(ax3, sc3, label='Pt Year')

        #Plot dz
        #
        c = dz
        vmin, vmax = malib.calcperc(c, perc=(5, 95))
        absmax = np.max(np.abs([vmin, vmax]))
        vmin = -absmax
        vmax = absmax
        im4 = ax4.imshow(hs_ma, cmap='gray', clim=hs_clim, alpha=0.5)
        sc4 = ax4.scatter(pX_fltr_mask_coreg, pY_fltr_mask_coreg, s=3, c=c, cmap='RdYlBu', vmin=vmin, vmax=vmax, edgecolors='none')
        cbar = pltlib.add_cbar(ax4, sc4, label='ICESat-GLAS - DEM (m)')

        for ax in (ax1, ax2, ax3, ax4):
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            ax.set_aspect('equal', 'box-forced')

        title='%s \n %i quality-filtered ICESat-GLAS shots, %i for valid surfaces, %i for co-registration' % ( pairname + ' (' + os.path.split(dem_fn)[1] +')', pX_fltr.shape[0], pX_fltr_mask_valsurf.shape[0], pX_fltr_mask_coreg.shape[0])
        fig.suptitle(title)
        fig.tight_layout()
        #This adjusts subplots to fit suptitle
        plt.subplots_adjust(top=0.92)
        #fig_fn = os.path.splitext(out_csv_fn)[0]+'.png'
        fig_fn = os.path.join(main_dir, pairname, pairname + "_" + os.path.splitext(os.path.split(out_csv_fn)[1])[0] +'.png')
        print "Saving figure: %s" % fig_fn
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight', pad_inches=0)
        plt.close(fig)