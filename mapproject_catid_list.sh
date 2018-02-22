#!/bin/bash
#
# Create Pan HRSI data that is QGIS-ready
#
# QUERY ADAPT dB by catid
# CORRECT & MOSAIC each P1BS scene to make a strip
# MAPPROJECT the P1BS strip using input DEM

t_start=$(date +%s)

function gettag() {
    xml=$1
    tag=$2
    echo $(grep "$tag" $xml | awk -F'[<>]' '{print $3}')
}

ADAPT=true
list_catid=$(tac ${1})
main_dir=$2
rpcdem=$3

#QUERY and return catid dirs with symlinks to data
if [ "$ADAPT" = true ] ; then

    echo; echo "Determine RPCDEM prj used to mapproject input prior to stereo ..."
    proj_rpcdem=$(proj_select.py ${rpcdem})

    for catid in $list_catid ; do
        mkdir -p ${main_dir}/${catid}
        cmd=''
        echo; echo "Querying ngadb, putting the symlinks catid ${catid} in ${main_dir}/${catid}"; echo
        cmd+="time query_db_catid.py $catid -out_dir ${main_dir}/${catid} ; "
        cmd_list+=\ \'$cmd\'
    done

    # Do the ADAPT db querying in parallel before any other processing
    #eval parallel --delay 2 -verbose -j 20 ::: $cmd_list

    # Loop over catid list again, do processing
    for catid in $list_catid ; do
        cd ${main_dir}/${catid}
        # ---FOR NOW -- only correct and mosaic P1BS data
        #ntfmos.sh ${main_dir}/${catid}

        mapproject_fill.py ${catid}.r100.tif $rpcdem

    done


else
    echo; echo "HTis workflow not running on ADAPT."; echo
fi



### Set up the filename you want
##mosMerge=${sensor}_${date}_${catID}.r100.M1BS.tif
##
### MOSAIC
### For each band of MS, mosaic all scenes associated with the strip of data (denoted with the catID)
### NOTE the '&' at the end of dg_mosaic, which launches the process in the background
##
##for bandNum in `seq 1 4`;
##do
##	dg_mosaic --band=$bandNum *${catID}_*M1BS*${ext} --output-prefix ${sensor}_${date}_${catID} --reduce-percent=100 &
##	# MAPPROJECT
##	# to Geographic WGS84 (EPSG:4326) -- you can choose another prj if you want.
##	mapproject --nodata-value=-99 --threads=20 -t rpc --t_srs 'EPSG:4326' $rpcdem $mosMerge ${mosMerge%.tif}.xml ${mosMerge%.tif}_prj.tif
##
##	# you are crazy if you don't build pyramid layers...
##	gdaladdo -ro -r average ${mosMerge%.tif}_prj.tif 2 4 8 16 32 64
##
##
##done
##
##wait
##
### STACK
### Stack the mosaiced bands
##gdal_merge.py -separate -o $mosMerge `ls ${sensor}_${date}_${catID}*tif`
##
### Now, COPY XML from a band's XML to the merged dataset
##cp ${sensor}_${date}_${catID}.r100.b1.xml ${mosMerge%.tif}.xml

t_end=$(date +%s)
t_diff=$(expr "$t_end" - "$t_start")
t_diff_hr=$(printf "%0.4f" $(echo "$t_diff/3600" | bc -l ))

echo; date
echo "Total processing time for pair ${pairname} in hrs: ${t_diff_hr}"