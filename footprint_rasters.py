#!/usr/bin/python

# Import and function definitions
import os, sys, osgeo, subprocess as subp
from osgeo import ogr, osr, gdal

import argparse
import dsm_info

def run_wait_os(cmdStr, print_stdOut=True):
    """
    Initialize OS command
    Wait for results (Communicate results i.e., make python wait until process is finished to proceed with next step)
    """
    import subprocess as subp

    Cmd = subp.Popen(cmdStr.rstrip('\n'), stdout=subp.PIPE, shell=True)
    stdOut, err = Cmd.communicate()

    if print_stdOut:
        print ("\tInitialized: %s" %(cmdStr))
        #print ("\t..Waiting for command to run...")
        print("\t" + str(stdOut) + str(err))
        print("\tEnd of command.")

def make_kml(fn):
    print("\tGenerating kml...")
    kml_fn = os.path.splitext(fn)[0]+'.kml'

    cmdStr = "ogr2ogr -f KML {} {}".format(kml_fn, fn)
    cmd = subp.Popen(cmdStr, stdout=subp.PIPE, shell=True)
    s,e = cmd.communicate()

def make_csv(fn):
    print("\tGenerating csv...")
    csv_fn = os.path.splitext(fn)[0]+'.csv'

    cmdStr = "ogr2ogr -f CSV {} {}".format(csv_fn, fn)
    cmd = subp.Popen(cmdStr, stdout=subp.PIPE, shell=True)
    s,e = cmd.communicate()

def make_link(pairname, rootDir):
    """Create a links in 3 top-level dirs to the 3 types of output (CLR, DRG, and DEM) tifs in native pairname dirs
    """

    src_tuple = (
                "out-DEM_24m_hs_az315.tif",\
				"out-DEM_4m_hs_az315.tif",\
				"out-DEM_1m_hs_az315.tif",\
                "_color_hs.tif",\
                "_ortho.tif",\
                "_ortho_4m.tif",\
                "out-DEM_1m.tif",\
                "out-DEM_4m.tif",\
                "out-DEM_24m.tif"\
                )

    # Find the file using the src string
    for root, dirs, files in os.walk(os.path.join(rootDir,pairname)):
        for fyle in files:
            if fyle.endswith(src_tuple):
                srcFull = os.path.join(rootDir,pairname,fyle)

                outDir = ''
                fyle = os.path.split(srcFull)[1]

                if not pairname in fyle and any (x in fyle for x in ["out-strip-DEM.tif", "out-DEM.tif", "out-DEM_1m.tif", "out-DEM_4m.tif", "out-DEM_24m.tif"]):
                    outDir = os.path.join(rootDir, "_dem")
                    dst = os.path.join(outDir, pairname+'_DEM'+fyle.split('-DEM')[1])

                if "color_hs.tif" in fyle:
                    outDir = os.path.join(rootDir, "_color_hs")
                    dst = os.path.join(outDir, pairname+fyle.split('DEM')[1])

                if "out-DEM_24m_hs_az315.tif" in fyle:
                    outDir = os.path.join(rootDir, "_hs","24m")
                    dst = os.path.join(outDir, pairname+fyle.split('DEM')[1])

                if "out-DEM_4m_hs_az315.tif" in fyle:
                    outDir = os.path.join(rootDir, "_hs","4m")
                    dst = os.path.join(outDir, pairname+fyle.split('DEM')[1])

                if "out-DEM_1m_hs_az315.tif" in fyle:
                    outDir = os.path.join(rootDir, "_hs","1m")
                    dst = os.path.join(outDir, pairname+fyle.split('DEM')[1])

                if pairname in fyle and "_ortho" in fyle:
                    outDir = os.path.join(rootDir, "_ortho")
                    dst = os.path.join(outDir, pairname+'_ortho'+fyle.split('_ortho')[1])

                if outDir:
                    os.system('mkdir -p %s' % outDir)

                    if os.path.isfile(dst):
                        os.remove(dst)

                    if os.path.isfile(srcFull):
                        cmdStr = "ln -s {} {}".format(srcFull, dst)
                        cmd = subp.Popen(cmdStr, stdout=subp.PIPE, shell=True)
                        print("\t\tWriting symlink: %s " %dst)


def getparser():
    parser = argparse.ArgumentParser(description="Create footprints of rasters files")
    parser.add_argument('ras_dir', default=None, help='Path to dir with raster to footprint')
    parser.add_argument('out_dir', default=None, type=str, help='Output dir out_shp')
    parser.add_argument('-ras_ext', default='.tif', help='The extension of rasters to be footprinted')
    parser.add_argument('-out_shp', default='raster_footprints', help='Output shapefile name of footprints')
    parser.add_argument('-file_fieldname', type=str, default='FILE', help='String indicating the field name describing the files')
    parser.add_argument('-c_pct', default='.25', type=str, help='The percent by which input pixel sizes will be coarsened (divided by)')
    parser.add_argument('-tmp_dir', default=None, type=str, help='Output dir for tmp files')
    parser.add_argument('-dir_exc_list', nargs='+', default=None, help='Exclude subdirs that start with strings in this list (_ ex z)')
    parser.add_argument('-dsm', action='store_true', default=False, help='footprint DSMs')
    parser.add_argument('-kml', action='store_true', default=False, help='Output kml of footprints for Google Earth')
    parser.add_argument('-csv', action='store_true', default=False, help='Output csv of attributes')
    parser.add_argument('-link', action='store_true', default=False, help='Write a symlink')
    parser.add_argument('-proj', type=str, default='polar', help='Specify prj type with a name')
    return parser

def main():

    """Creates/Update a top-level directory's footprint shapefile of all rasters meeting a specified extension.

        Produces a coarsened shapefile of the valid pixels from all these rasters with file and path attributes.

        The -dsm flag returns a shapefile that includes many attributes (related to the stereo acquisition) from the associated XMLs found in the same dir as the raster.
        The -kml and -csv flags will return a KML or CSV version of the shapefile

    Note: A particular version of SQLite is needed to run the SQL Gunion operation
    This version should be sourced prior to running this script, by appending yout PATH variable to dirs that hold the correct version of SQLite.
    I have this:
        export PATH=/usr/local/bin:/usr/bin:/bin:/opt/StereoPipeline/bin:/opt/StereoPipeline/libexec:/opt/bin:/opt/exelis/idl/bin:/usr/mpi/gcc/openmpi-1.10.5a1/bin:/opt/PGSC-imagery_utils:$PATH
    in this file:
        $HOME/code/sqlite_fix_new.env
    and source like this:
        source $HOME/code/sqlite_fix_new.env
    """
    parser = getparser()
    args = parser.parse_args()

    ras_dir = args.ras_dir
    out_dir = args.out_dir
    ras_ext = args.ras_ext
    out_shp = args.out_shp
    file_fieldname = args.file_fieldname
    c_pct = args.c_pct
    tmp_dir = args.tmp_dir
    dir_exc_list = args.dir_exc_list
    DSM = args.dsm
    KML = args.kml
    CSV = args.csv
    LINK = args.link

    if tmp_dir is None:
        tmp_dir = out_dir

    out_shp_fn = os.path.join(out_dir,out_shp)

    if not out_shp.endswith('shp'):
        out_shp_fn += '.shp'

    # Projections
    proj = args.proj
    #default is polar
    proj4 = "'+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'"
    
    if 'modis' in proj:
        proj4 = "'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'"
        proj4 = "'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'"
    if 'wgs' in proj:
        proj4 = "'+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'"

    print "\n\tRunning footprints on: %s\n" %ras_dir

    # Collect raster feature names in working directory and subfolders therein
    ras_name_list = []
    pathroot = []
    ras_fn_list = []

    print "\tWalking main dir, building list of rasters..."

    for root, dirs, files in os.walk(ras_dir):
        ## https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
        if dir_exc_list is not None:
            dirs[:] = [d for d in dirs if not d.startswith(tuple(dir_exc_list))]
        for f in files:
            if f.endswith(ras_ext) or f.endswith(ras_ext.upper()):
                ras_fn_list.append(os.path.join(root, f))
                ras_name_list.append(f)
                pathroot.append(root)

    print "\tIterating over raster list..."
    for num, ras_fn in enumerate(ras_fn_list):

        # Clean up tmp files which may otherwise interfere with footprinting
        file_list = os.listdir(tmp_dir)
        for f in file_list:
            if 'tmp' in f:
                os.remove(os.path.join(tmp_dir,f))

        path_name,file_name = os.path.split(ras_fn)
        dir_name = os.path.split(path_name)[1].replace('.','_').replace('-','_')
        tmp_file_name = dir_name + "_" + file_name.strip('.tif').replace('.','_').replace('-','_')
        tmp1 = os.path.join(tmp_dir, "tmp1_"+tmp_file_name+".tif")
        tmp2 = os.path.join(tmp_dir, "tmp2_"+tmp_file_name+".tif")
        tmp3 = os.path.join(tmp_dir, "tmp3_"+tmp_file_name+".shp")
        tmp4 = os.path.join(tmp_dir, "tmp4_"+tmp_file_name+".shp")
        tmp5 = os.path.join(tmp_dir, "tmp5_"+tmp_file_name+".shp")
        tmp6 = os.path.join(tmp_dir, "tmp6_"+tmp_file_name+".shp")
        tmp_final = os.path.join(tmp_dir, "tmp_final_"+tmp_file_name+".shp")

        if not os.path.isfile(ras_fn):
            print "\tWill not footprint: %s does not exist." %ras_fn
        else:

            print "\n\t # %s of %s rasters " %(num+1,len(ras_fn_list))
            print "\tOutput shp: %s" %out_shp_fn
            try:
                print "\tChecking: %s" %ras_fn
                # Check to see if the file name exists in the 'Name' field of the output shp
                update = True
                path_fieldname = 'PATH'
                if os.path.isfile(out_shp_fn):
                    # Open out_shp_fn and get field names"
                    shp = ogr.Open(out_shp_fn, 1)
                    layer = shp.GetLayer()
                    layer_defn = layer.GetLayerDefn()
                    field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]

                    for feature in layer:
                        if feature.GetField(file_fieldname) == file_name and feature.GetField(path_fieldname) == path_name:
                            print "\tFootprint of %s already exists" %file_name
                            update = False
                            break
                        pass
                    shp, layer, layer_defn, field_names, feature = (None for i in range(5))

                if update:

                    print "\tFootprinting: %s" %os.path.split(ras_fn)[1]
                    print "\t\t...COARSEN..."
                    cmdStr = "gdal_translate -outsize {}% {}% -co compress=lzw -b 1 -ot Float32 {} {}".format(c_pct, c_pct, ras_fn, tmp1)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...CALC BINARY MASK..."
                    cmdStr = 'gdal_calc.py --overwrite -A {} --outfile={} --calc="1*((A>=-50)+0*(A<-50))" --NoDataValue=-99'.format(tmp1, tmp2)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...SIEVE..."
                    cmdStr = 'gdal_sieve.py -8 -st 5 {} {}'.format(tmp2, tmp3)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...POLYGONIZE..."
                    cmdStr = "gdal_polygonize.py -8 -f 'ESRI Shapefile' {} {}".format(tmp3, tmp4)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...REMOVE NODATA POLYGONS..."
                    cmdStr = "ogr2ogr {} {} -where 'DN>0'".format(tmp5, tmp4)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...REPROJECT..."
                    cmdStr = "ogr2ogr -f 'ESRI Shapefile' -t_srs {} {} {} -overwrite".format(proj4,tmp6, tmp5)
                    run_wait_os(cmdStr,print_stdOut=False)
                    print "\t\t...DISSOLVE/AGGREGATE INTO 1 FEATURE..."
                    input_basename = os.path.split(tmp6)[1].replace(".shp","")
                    cmdStr = "ogr2ogr {} {} -dialect sqlite -sql 'SELECT GUnion(geometry), DN FROM {} GROUP BY DN'".format(tmp_final, tmp6, input_basename)
                    run_wait_os(cmdStr,print_stdOut=False)

                if update and os.path.isfile(tmp_final):
                    # Add fields to shp
                    ##https://gis.stackexchange.com/questions/3623/how-to-add-custom-feature-attributes-to-shapefile-using-python
                    # Open a Shapefile, and get the layer
                    shp = ogr.Open(tmp_final, 1)
                    layer = shp.GetLayer()

                    # [1] Set the lists for field names and attributes
                    field_names_list = [file_fieldname,path_fieldname]
                    field_attributes_list = [os.path.split(ras_fn)[1], os.path.split(ras_fn)[0]]

                    if DSM:
                        #print "Adding 'pairname' field and attribute before getting DSM Info."
                        pairname_fieldname = ['PAIRNAME']
                        field_names_list += pairname_fieldname
                        pairname = os.path.split(os.path.split(ras_fn)[0])[1]
                        field_attributes_list += [pairname]

                        # Kick out ang_conv (c), ang_bie (b), ang_asm (a), the DSM info header comma-delim'd string (dsm_hdr), and the attributes comm-delim'd string associated with that header
                        c,b,a,dsm_hdr,attributes = dsm_info.main(path_name)

                        #print "\t\tHeader: %s" %(dsm_hdr)
                        #print "\t\tAttributes: %s" %(attributes)

                        dsm_hdr_list = dsm_hdr.rstrip().strip(',').split(',')
                        attributes_list = attributes.rstrip().strip(',').split(',')

                        field_names_list += dsm_hdr_list
                        field_attributes_list += attributes_list

                    # [2] Enumerate field_names_list to get corresponding attribute from field_attributes_list, then get each attribute type, and set type of each new field
                    for num, new_field_name in enumerate(field_names_list):
                        if any(x in field_attributes_list[num] for x in ['_','/',':']):
                            fieldType = ogr.OFTString
                        elif any(x in field_names_list[num] for x in ['YEAR','MONTH','DOY']):
                            fieldType = ogr.OFTInteger
                        else:
                            fieldType = ogr.OFTReal

                        new_field = ogr.FieldDefn(new_field_name, fieldType)
                        layer.CreateField(new_field)

                    # [3] From field_attributes_list, update fields with attributes
                    ## http://www.digital-geography.com/create-and-edit-shapefiles-with-python-only/#.V8hlLfkrLRY
                    layer_defn = layer.GetLayerDefn()
                    field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
                    feature = layer.GetFeature(0)   # Gets first, and only, feature
                    for feature in layer:
                        for num, f_name in enumerate(field_names):
                            if num > 0: ## the 0 idx is the useless 'DN' field
                                feature.SetField(num, field_attributes_list[num-1])
                                layer.SetFeature(feature)
                    shp = None

                    # Append final tmp to out_shp
                    if os.path.isfile(out_shp_fn):
                        print "\tUpdating footprint: %s" %out_shp_fn
                        cmdStr = "ogr2ogr -f 'ESRI Shapefile' -update -append {} {}".format(out_shp_fn, tmp_final)
                        run_wait_os(cmdStr,print_stdOut=False)
                    else:
                        print "\tCreating footprint: %s" %out_shp_fn
                        cmdStr = "ogr2ogr -f 'ESRI Shapefile' {} {}".format(out_shp_fn, tmp_final)
                        run_wait_os(cmdStr,print_stdOut=False)

                    if DSM and LINK:
                        print "\tCreating symlinks: %s" %pairname
                        make_link(pairname, os.path.split(os.path.split(ras_fn)[0])[0])

                # Clean up tmp files
                file_list = os.listdir(tmp_dir)
                for f in file_list:
                    if 'tmp' in f:
                        os.remove(os.path.join(tmp_dir,f))

            except Exception, e:
                print "\tFailed to footprint (are pairname XML files present?) : %s" %ras_fn

    if KML and os.path.isfile(out_shp_fn):
        make_kml(out_shp_fn)
    if CSV and os.path.isfile(out_shp_fn):
        make_csv(out_shp_fn)

if __name__ == '__main__':
    main()