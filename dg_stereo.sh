#!/bin/bash
#
# DEM Workflow: wv_correct, mosaic, mapproject, stereo, point2dem, hillshades, & orthoimages for individual stereopairs on DISCOVER & ADAPT
# paul montesano, david shean, maggie wooten, christopher neigh
#
# pairname=WV01_20140603_102001002F42A400_1020010031EBEF00
# example of call on DISCOVER:
#     dg_stereo.sh $pairname false
# example of call on ADAPT:
#	  pupsh "hostname ~ 'himat115'" "dg_stereo.sh WV02_20160512_10300100548DD500_103001005422BA00 true true false true batch_andes '' true $HOME/code/nodes_ecotone_ngaproc false 21 1024"
#	    or
#     dg_stereo.sh $pairname true
#       or
#     pupsh "hostname ~ 'ecotone16'" "dg_stereo_par.sh /att/pubrepo/DEM/hrsi_dsm/list_pairname"
#
# Dependencies (sh & python scripts run as cmd line tools):
#   query_db_catid.py		script that returns the ADAPT dir of images of given catid
#   proj_select.py          get the best prj used to mapproject input
#   utm_proj_select.py		force get UTM prj for DEM and ortho; edit to script from pygeotools to force select UTM zone (instead of best prj); 
#   color_hs.py
#   ntfmos.sh
#   dg_stereo_int.py
#   warptool.py
# Note: on ADAPT, will run parallel_stereo on launch node only, thus, no nodeslist needed.


t_start=$(date +%s)

function gettag() {
    xml=$1
    tag=$2
    echo $(grep "$tag" $xml | awk -F'[<>]' '{print $3}')
}

#Hardcoded Args
tile_size=1024 #20480

# Required Args
pairname=$1
TEST=$2           #true or false
ADAPT=$3          #true or false
MAP=$4            #true or false
RUN_PSTEREO=$5    #true or false
batch_name=$6
rpcdem=$7         #can be blank var ''
NODES=$8          #true or false
nodeslist=$9
SGM=${10}         #true or false
# Using 7 for vegetation (very noisy, but resolves more gaps?) - use 21+ for mountains/glaciers/other terrain
subpix_kern=${11}
# 1024 probably decent
erode_max_size=${12}
# 21 is default
corr_kern=${13}
# 300 is default; increase for more difficult areas.
corr_time=${14} 

if [ "$ADAPT" = false ]; then
    TEST=false
fi

if [ "$TEST" = true ]; then
    # Optional Args (stereogrammetry testing)
    crop=${15}   #"0 190000 40000 40000"
    #sa=${16}	   #if sgm is true, then use 1 for sgm or 2 for mgm
    #cm=${17}      #cost mode for stereo
fi

if [ "$ADAPT" = true ]; then
    out_root=$NOBACKUP/outASP/${batch_name}
    mkdir -p $out_root
else
    #out_root=$4 # output directory is 4th input if on DISCOVER
    out_root=/discover/nobackup/projects/boreal_nga/ASP/${batch_name}
fi

left_catid="$(echo $pairname | awk -F '_' '{print $3}')"
right_catid="$(echo $pairname | awk -F '_' '{print $4}')"

if [ -z "$SGM" ]; then
	SGM=false
fi

ncpu=$(cat /proc/cpuinfo | egrep "core id|physical id" | tr -d "\n" | sed s/physical/\\nphysical/g | grep -v ^$ | sort | uniq | wc -l)

gdal_opts="-co TILED=YES -co COMPRESS=LZW -co BIGTIFF=YES"
gdal_opts+=" -co BLOCKXSIZE=256 -co BLOCKYSIZE=256"
#gdal_opts+=" -co NUM_THREADS=$ncpu"

parallel_point2dem=false
if [ "$ncpu" -gt "12" ] ; then
    parallel_point2dem=true
fi

# Stereogrammetry
out=${out_root}/${pairname}/out
stereo_opts=''
stereo_args=''
sgm_opts=''

if [ -e ${out}-strip-PC.tif ]; then
	mv ${out}-strip-PC.tif ${out}-PC.tif
fi
if [ -e ${out}-strip-DEM.tif ]; then
	mv ${out}-strip-DEM.tif ${out}-DEM_native.tif
fi

#Set entry point based on contents of outdir
if [ -e ${out}-PC.tif ]; then
    e=5
elif [ -e ${out}-F.tif ]; then
    e=4
elif [ -e ${out}-RD.tif ]; then
    e=3
elif [ -e ${out}-D.tif ]; then
    e=2
elif [ -e ${out}-R_sub.tif ]; then
    e=1
else
    e=0
fi

#Set in_left and in_right consistent with expected output of dg_mosaic
in_left=${out_root}/${pairname}/${left_catid}.r100.tif
in_right=${out_root}/${pairname}/${right_catid}.r100.tif

#Set the name of the output orthoimage
ortho_ext=_ortho.tif
out_ortho=${out_root}/${pairname}/${pairname}${ortho_ext}

if [ ! -e $in_left ] || [ ! -e $in_right  ] ; then
    mkdir -p ${out_root}/${pairname}
    if [ ! -e ${out_ortho} ] ; then
        if [ "$ADAPT" = true ] ; then
            for catid in $left_catid $right_catid ; do
                cmd=''
                echo; echo "Querying ngadb, putting the symlinks catid ${catid} in ${out_root}/${pairname}"; echo
                cmd+="time query_db_catid.py $catid -out_dir ${out_root}/${pairname} ; "
                cmd_list+=\ \'$cmd\'
            done

            # Do the ADAPT db querying in parallel
        
            eval parallel --delay 2 -verbose -j 2 ::: $cmd_list
        else
            echo; echo "Workflow not running on ADAPT, querying for input already done."; echo
        fi
    fi
fi

if [ ! -e "${out}-PC.tif" ] ; then
    echo; echo "Running wv_correct and dg_mosaic to create:"; echo "${in_left}"; echo "${in_right}"
    ntfmos.sh ${out_root}/${pairname}
fi

if [ ! -e $in_left ] && [ ! -e ${in_left%.*}.xml ]; then
    in_left_xml=$(echo $(ls ${out_root}/${pairname}/*${left_catid}*P1BS*.xml | grep -v aux | head -1))
else
    in_left_xml=${in_left%.*}.xml
fi
if [ ! -e $in_right ] && [ ! -e ${in_right%.*}.xml ]; then
    in_right_xml=$(echo $(ls ${out_root}/${pairname}/*${right_catid}*P1BS*.xml | grep -v aux | head -1))
else
    in_right_xml=${in_right%.*}.xml
fi

# Get proj from XML
if [ "$MAP" = true ] ; then
    echo; echo "Determine RPCDEM prj used to mapproject input prior to stereo ..."
    proj_rpcdem=$(proj_select.py ${rpcdem})
fi

echo; echo "Determine output UTM prj, and native resolution ..."
proj=$(utm_proj_select.py ${in_left_xml})
echo "Projection: ${proj}"

if grep -q MEANPRODUCTGSD $in_left_xml ; then
    res1=$(printf '%.3f' $(gettag $in_left_xml 'MEANPRODUCTGSD'))
else
    res1=$(printf '%.3f' $(gettag $in_left_xml 'MEANCOLLECTEDGSD'))
fi
if grep -q MEANPRODUCTGSD $in_right_xml ; then
    res2=$(printf '%.3f' $(gettag $in_right_xml 'MEANPRODUCTGSD'))
else
    res2=$(printf '%.3f' $(gettag $in_right_xml 'MEANCOLLECTEDGSD'))
fi
echo "GSD resolutions"
echo "${left_catid}: $res1 GSD"
echo "${right_catid}: $res2 GSD"

if [ $(echo "a=($res1 < $res2); a" | bc -l) -eq 1 ] ; then
    native_res=$res1
    echo "Native res is from $left_catid : ${native_res}"
    mos4ortho_img=$in_left
    mos4ortho_catid=$left_catid
else
    native_res=$res2
    echo "Native res is from $right_catid : ${native_res}"
    mos4ortho_img=$in_right
    mos4ortho_catid=$right_catid
fi

if [ "$e" -lt "5" ] && [ -e $in_left ] && [ -e $in_right ] ; then
    stereo_opts+="-t dg"

    #Map mosaiced input images using ASP mapproject
    if [ "$MAP" = true ] ; then
        map_opts="--threads $ncpu -t rpc --nodata-value 0 --t_srs \"$proj_rpcdem\""

        if [[ -n $native_res ]]; then
            map_opts+=" --tr $native_res"
            outext="${outext}_${native_res}m"
        fi
        echo "Projection used for initial alignment of stereopairs:"
        echo $proj_rpcdem
        echo "Computing intersection extent in projected coordinates:"
        map_extent=$(dg_stereo_int.py $in_left_xml $in_right_xml "$proj_rpcdem")
        echo $map_extent; echo

        for in_img in $in_left $in_right; do
            ln -sv ${in_img%.tif}.xml ${in_img%.tif}${outext}.xml
            map_arg="--t_projwin $map_extent $rpcdem ${in_img} ${in_img%.tif}${outext}.xml ${in_img%.tif}${outext}.tif"
            if [ ! -e ${in_img%.tif}${outext}.tif ]; then
                echo; date; echo;
                echo mapproject $map_opts $map_arg
                eval time mapproject $map_opts $map_arg
            fi
        done
        rpcdem_warp=${out_root}/${pairname}/$(basename ${rpcdem%.*})_warp.tif
        if [ ! -e $rpcdem_warp ] ; then
            echo; echo "Clip the VRT rpcdem with the mapprojected extent..."; echo
            warptool.py -tr 'last' -te 'first' ${in_img%.tif}${outext}.tif $rpcdem -outdir ${out_root}/${pairname}
        fi
        # Rename rpcdem to the clipped file
        rpcdem=$rpcdem_warp

        stereo_args+="$rpcdem"
        stereo_opts+=" --alignment-method None"

    #Don't map inputs, let ASP do the alignment
    else
        echo; date; echo;
        outext=""
        stereo_opts+=" --alignment-method AffineEpipolar"
    fi

    stereo_opts+=" --subpixel-kernel $subpix_kern $subpix_kern"
    stereo_opts+=" --erode-max-size $erode_max_size"
    stereo_opts+=" --corr-kernel $corr_kern $corr_kern"
    stereo_opts+=" --corr-timeout $corr_time"
    stereo_opts+=" --individually-normalize"
    stereo_opts+=" --tif-compress LZW"

    if [ ! -z "$crop" ]; then
        stereo_opts+=" --left-image-crop-win $crop"
    fi

    # Done like this so, if present, rpcdem is last
    stereo_args="${in_left%.*}${outext}.tif ${in_right%.*}${outext}.tif ${in_left%.*}${outext}.xml ${in_right%.*}${outext}.xml ${out} $stereo_args"

    # ADAPT processing needs these
    par_opts="--job-size-w $tile_size --job-size-h $tile_size"
    par_opts+=" --threads-singleprocess $ncpu"
    par_opts+=" --processes $ncpu"
    par_opts+=" --threads-multiprocess 1"
    if [ "$NODES" = true  ] ; then
        par_opts+=" --nodes-list=$nodeslist"
    fi

    if [ "$RUN_PSTEREO" = true ] && [ "$SGM" = true ] ; then
        # SGM stereo runs. Not applicable for our DISCOVER processing
        if [ ! -z "$sa" ]; then
            sgm_opts+=" --stereo-algorithm $sa"
        else
            sgm_opts+=" --stereo-algorithm 1"
        fi
        if [ ! -z "$cm" ]; then
            sgm_opts+=" --cost-mode $cm"
        else
            sgm_opts+=" --cost-mode 3"
        fi
        sgm_opts+=" --corr-memory-limit-mb 5000"
        sgm_opts+=" --corr-kernel 7 7"
        sgm_opts+=" --corr-tile-size $tile_size"
        sgm_opts+=" --xcorr-threshold -1"
        sgm_opts+=" --subpixel-mode 0"
        sgm_opts+=" --median-filter-size 3"
        sgm_opts+=" --texture-smooth-size 7"
        sgm_opts+=" --texture-smooth-scale 0"
        sgm_opts+=" --threads 10"
        sgm_opts+=" $stereo_opts"

        echo; date; echo;
        #eval time stereo -e $e $sgm_opts $stereo_args
        cmd="parallel_stereo -e $e $par_opts $sgm_opts $stereo_args"
        echo $cmd
        eval $cmd
    else
        # DISCOVER processing needs these.
        #stereo_opts+=" --corr-kernel $corr_kern $corr_kern" #mw 7/2 this was causing error in stereo call. (calling --corr-kernel twice in the command). temporarily comment
        stereo_opts+=" --subpixel-mode 2"
        stereo_opts+=" --filter-mode 1"
        stereo_opts+=" --cost-mode 2"

        if [ "$RUN_PSTEREO" = true ] && [ "$ADAPT" = true ] ; then
            echo; echo $stereo_args ; echo
            echo; echo "parallel_stereo -e $e $par_opts $stereo_opts $stereo_args"; echo
            eval time parallel_stereo -e $e $par_opts $stereo_opts $stereo_args
            echo; echo "Removing intermediate logs..."
            rm ${out}-log-stereo_parse*.txt
        else
            echo; echo "stereo $stereo_opts $stereo_args"; echo
            eval time stereo -e $e $stereo_opts $stereo_args
        fi
    fi
fi
if [ ! -e "${out}-PC.tif" ] ; then
    echo; echo "Stereogrammetry unsuccessful. Exiting."
    exit 1
else
    echo; echo "Stereo point-cloud file exists."
    if [ "$ADAPT" = true ] && gdalinfo ${out}-PC.tif | grep -q VRT ; then
        echo; echo "Convert PC.tif from virtual to real"; echo
        eval time gdal_translate $gdal_opts ${out}-PC.tif ${out}-PC_full.tif
        mv ${out}-PC_full.tif ${out}-PC.tif
    fi
    if [ "$RUN_PSTEREO" = true ] ; then
        echo; echo "Removing intermediate parallel_stereo dirs..."
        echo "${out}-"*/
        rm -rf ${out}-*/
    fi

    stats_res=24
    mid_res=4
    fine_res=1

    stats_dem=${out}-DEM_${stats_res}m.tif
    mid_dem=${out}-DEM_${mid_res}m.tif
    fine_dem=${out}-DEM_${fine_res}m.tif

    # DEM Generation
    cmd_list=''
    dem_ndv=-99

    base_dem_opts=" --remove-outliers --remove-outliers-params 75.0 3.0"
    base_dem_opts+=" --threads 4"
    base_dem_opts+=" --t_srs \"$proj\""
    echo; echo "Check for dems..."; echo
    for dem_res in $stats_res $mid_res $fine_res ; do
        dem_opts="$base_dem_opts"
        if [ ! -e ${out}-DEM_${dem_res}m.tif ]; then
            echo "Creating DEM at ${dem_res}m ..."
            dem_opts+=" --nodata-value $dem_ndv"
    	    dem_opts+=" --tr $dem_res"
            if [ "$dem_res" = "$stats_res" ] ; then
                dem_opts+=" --dem-hole-fill-len 10"
            fi
            dem_opts+=" -o ${out}_${dem_res}m"

          echo; date; echo;
          echo point2dem $dem_opts ${out}-PC.tif
          echo

          if [ "$parallel_point2dem" = true ] ; then
              cmd=''
              cmd+="time point2dem $dem_opts ${out}-PC.tif; "
              cmd+="mv ${out}_${dem_res}m-DEM.tif ${out}-DEM_${dem_res}m.tif; "
              cmd_list+=\ \'$cmd\'
          else
              eval time point2dem $dem_opts ${out}-PC.tif
              mv ${out}_${dem_res}m-DEM.tif ${out}-DEM_${dem_res}m.tif
          fi
        else
            echo "Finished: ${out}-DEM_${dem_res}m.tif"
        fi
        
    done

    if [[ ! -z $cmd_list ]] ; then

       	if (( $ncpu > 15 )) ; then
            njobs=4
        else
            njobs=2
        fi
        eval parallel -verbose -j $njobs ::: $cmd_list
    fi

    # Color Shaded Relief Generation
    mean=$(gdalinfo -stats $stats_dem | grep MEAN | awk -F '=' '{print $2}')
    stddev=$(gdalinfo -stats $stats_dem | grep STDDEV | awk -F '=' '{print $2}')

    min=$(echo $mean $stddev | awk '{print $1 - $2}')
    max=$(echo $mean $stddev | awk '{print $1 + $2}')

    cmd_list=''
    echo; echo "Check for color-shaded reliefs..."; echo
    for dem in $stats_dem $mid_dem $fine_dem ; do
    	cmd=''
        if [ ! -e ${dem%.*}_color_hs.tif.ovr ]; then
            rm -f ${dem%.*}_color_hs.tif

    	    cmd+="time color_hs.py $dem -clim $min $max -hs_overlay -alpha .8; "
    	    cmd_list+=\ \'$cmd\'
        else
            echo "Finished: ${dem%.*}_color_hs.tif"
        fi
    done
    if [[ ! -z $cmd_list ]]; then
        echo; echo "Do all colorshades and hillshades in parallel"; echo
        eval parallel -verbose -j 10 ::: $cmd_list
        rm -v ${out}*color.tif
    fi

    ortho_opts="--nodata-value 0"

    map_opts="$ortho_opts"
    map_opts+=" -t rpc"
    map_opts+=" --num-processes $ncpu"

    #If both in_left and in_right exist, then catid mosaics are complete, and in_left can be ortho'd
    # else no mosiacs done, in_left is an xml used for proj and native_res; need indiv scenes indiv ortho'd then dem_mosaic
    if [ ! -e ${out_ortho} ] ; then
        if [ -e ${mos4ortho_img} ] ; then
            echo; echo "Mapproject at ${native_res}m ${mos4ortho_img} onto ${stats_dem}"; echo
            map_opts=" --tr $native_res"
            map_args="$stats_dem $mos4ortho_img ${mos4ortho_img%.*}.xml ${out_ortho}"
            time mapproject $map_opts $map_args
        else
            # This case exists to handle pairname dirs that dont have *.r100.tif; so, for each ntf run mapprj then use dem_mosaic
            echo; echo "Mapproject each indiv NTF onto ${stats_dem}"; echo
            echo "CATID: ${mos4ortho_catid}"
            ntf_list=$(ls ${out_root}/${pairname} | grep -e "${mos4ortho_catid}" | grep -i P1BS | egrep 'ntf|tif' | grep -v 'corr')
            
            if [ ! "$ntf_list" ] && [ "$ADAPT" = true ]  ; then
                 echo; echo "Get ADAPT dir with imagery to mapproject"; echo
                 query_db_catid.py ${mos4ortho_catid} -out_dir ${out_root}/${pairname}
                 ntf_list=$(ls ${out_root}/${pairname} | grep -e "${mos4ortho_catid}" | grep -i P1BS | egrep 'ntf|tif' | grep -v 'corr')
            fi

            cmd_list=''

            for ntf in $ntf_list ; do
                ntf_fn=${out_root}/${pairname}/${ntf}
                indiv_ortho=${out_root}/${pairname}/${ntf%.*}${ortho_ext}
                map_args="$stats_dem ${ntf_fn} ${ntf_fn%.*}.xml $indiv_ortho"
    	        echo $ntf_fn
    	        cmd=''
    	        cmd+="time mapproject $map_opts $map_args; "
                cmd_list+=\ \'$cmd\'
            done

            echo; echo "Do orthos for each P1BS scene running mapproject in parallel"; echo
            eval parallel -verbose -j 6 ::: $cmd_list

            echo; echo "Do dem_mosaic at native res of orthos"; echo
            echo; echo "dem_mosaic --tr $native_res --threads $ncpu `ls ${out_root}/${pairname}/*${ortho_ext}` -o ${out_root}/${pairname}/${pairname}"; echo

            time dem_mosaic --tr $native_res --threads $ncpu `ls ${out_root}/${pairname}/*${ortho_ext}` -o ${out_root}/${pairname}/${pairname}
            mv ${out_root}/${pairname}/${pairname}-tile-0.tif ${out_ortho}
        fi
    fi
    if [ ! -e ${out_ortho%.*}.tif.ovr ] ; then
        echo "Building overviews in background for:"; echo "$out_ortho"
        gdaladdo -ro -r average ${out_ortho} 2 4 8 16 32 64 &
    fi
    if [ ! -e ${out_ortho%.*}_${mid_res}m.tif.ovr ] ; then
        echo; echo "Simple gdal_translate to coarsen native_res ortho..."; echo
        gdal_translate -tr ${mid_res} ${mid_res} ${out_ortho} ${out_ortho%.*}_${mid_res}m.tif
        gdaladdo -ro -r average ${out_ortho%.*}_${mid_res}m.tif 2 4 8 16 32 64 &
    fi
    if [ "$ADAPT" = true ] ; then
        echo; echo "Writing symlinks"; echo
        for i in $out_ortho ${out_ortho%.*}_${mid_res}m.tif ; do
            ln -sfv ${i} ${out_root}/_ortho/$(basename ${i})
        done
        for i in $stats_dem $mid_dem $fine_dem ; do
            dembase=$(basename ${i})
            ln -sfv ${i} ${out_root}/_dem/${pairname}_${dembase:4}

            if [ -e ${dembase%.*}_color_hs.tif ] ; then
                colorbase=${dembase%.*}_color_hs.tif
                ln -sfv ${i%.*}_color_hs.tif ${out_root}/_color_hs/${pairname}_${colorbase:4}
            fi
        done
    fi

    if [ "$TEST" = true ] ; then
        echo; echo "TEST run: Keeping all intermediate files..."
    else
        echo; echo "Removing intermediate files..."
    
        if [ -e ${out}-DEM_native.tif ]; then
	 	    rm -v ${out}-DEM_native.tif
        fi

        for i in $(ls ${out_root}/${pairname}/*P1BS*ortho.tif); do
            rm -v ${i%.*}*
        done
        #for i in $(ls ${out_root}/${pairname}/*_corr.*); do
        rm -rf $(ls ${out_root}/${pairname}/*.{tif,xml} | grep _corr)
        rm ${out}-log-stereo_parse*.txt

        if [ -e "${out_ortho}" ] ; then
            rm -v "${out_root}/${pairname}/"*.r100*.tif
        fi

        rm -v "${out_root}/${pairname}/"out.*
        rm -v "${out_root}/${pairname}/"*warp.tif

        for i in sub.tif Mask.tif .match .exr center.txt ramp.txt; do
            rm -v ${out}-*${i}
        done
        for i in F L R RD D GoodPixelMap DEM-clr-shd DEM-hlshd-e25 DRG; do
            if [ -e "${out}-${i}.tif" ]; then
                rm -v "${out}-${i}.tif"
            fi
            if [ -e "${out}-strip-${i}.tif" ]; then
                rm -v "${out}-strip-${i}".*
            fi
        done
    fi
fi

t_end=$(date +%s)
t_diff=$(expr "$t_end" - "$t_start")
t_diff_hr=$(printf "%0.4f" $(echo "$t_diff/3600" | bc -l ))

echo; date
echo "Total processing time for pair ${pairname} in hrs: ${t_diff_hr}"
exit 1
