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
        if [ ! -e ${main_dir}/${catid}/${catid}_ortho_rpc.tif ] ; then
            
            echo "Running ntfmos.sh to apply wv_correct to scenes, and the mosaic scenes into strips for input into mapproject."
            ntfmos.sh ${main_dir}/${catid}
            
            echo "Running mapproject to ortho input strip."
            cmd=$(mapproject_fill.py ${main_dir}/${catid}/${catid}.r100.tif $rpcdem)
            echo $cmd
            eval $cmd
        fi
        gdaladdo -ro -r average ${main_dir}/${catid}/${catid}_ortho_rpc.tif 2 4 8 16 32 64
    done


else
    echo; echo "This workflow not running on ADAPT."; echo
fi

t_end=$(date +%s)
t_diff=$(expr "$t_end" - "$t_start")
t_diff_hr=$(printf "%0.4f" $(echo "$t_diff/3600" | bc -l ))

echo; date
echo "Total processing time for pair ${pairname} in hrs: ${t_diff_hr}"