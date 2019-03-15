#!/bin/bash
#
# Co-register a single HRSI DSM to ICESat-GLAS
#
# 3 steps
# - Masking input DEM to get control surfaces (dem_control.py)
# - Filtering GLAS over the control surfaces
# - PC aligning input DEM with filtered GLAS

source ~/anaconda3/bin/activate py2 #For Paul
#source ~/miniconda2/bin/activate py2 #My code

hostN=`/bin/hostname -s`
script_name=$(basename ${0})
pairname=${1}
main_dir=${2:-''}
dem_file=${3:-'out-DEM_1m.tif'}
glas=${4:-'/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/glas/circ_boreal/gla01-boreal50up-data.csv'}
ref_dem=${5:-'/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/tandemx/TDM90/mos/TDM1_90m_circ_DEM.vrt'}

elev_range="-15 15"
test_dir='/att/gpfsfs/briskfs01/ppl/pmontesa/tmp/test2'

if [ -z $main_dir ] ; then
    main_dir=$test_dir
fi

workdir=$main_dir/$pairname
mkdir -p $workdir
dem=$workdir/${dem_file}

#Pubrepo files needed: out-DEM_*.tif, *.xml; Need to needed files to NOBACKUP
ln -sf /att/pubrepo/DEM/hrsi_dsm/$pairname/*ortho*tif $main_dir/$pairname
xml_fn_list=$(ls /att/pubrepo/DEM/hrsi_dsm/$pairname/*.xml)
ln -sf $xml_fn_list $main_dir/$pairname
ln -sf /att/pubrepo/DEM/hrsi_dsm/$pairname/out-DEM*m.tif $main_dir/$pairname
if [ $main_dir = $test_dir ] ; then
    ln -sf /att/pubrepo/DEM/hrsi_dsm/$pairname/out-DEM*m*hs*.tif $main_dir/$pairname
fi

mkdir -p $main_dir/logs_coreg2glas
logfile=$main_dir/logs_coreg2glas/${script_name%.*}_${hostN}_${pairname}.log

echo; echo "Script call:" | tee $logfile
echo "${script_name} ${1} ${2} ${3}"| tee -a $logfile
echo "Main dir: $main_dir" | tee -a $logfile

dem_4control=${dem%1*}4m.tif
echo $dem_4control

echo; echo "Get static control DEM: masks input using TOA (dark & smooth) & DEM (roughness & slope) masked out of the 4m version of each DEM" ; echo | tee -a $logfile
cmd="dem_control.py $dem_4control -filt_param $ref_dem $elev_range"
echo $cmd ; echo  | tee -a $logfile
eval $cmd  | tee -a $logfile

echo; echo "Filter ICESat-GLAS (GLA14) using the *control.tif: this provides a set of GLAS over control surfaces only."
#using dem_4control, which expects a corresponding *_control.tif from dem_control.py
cmd="filter_glas_control.py $glas $dem_4control"
echo $cmd ; echo  | tee -a $logfile
eval $cmd  | tee -a $logfile

glas_file=$(basename $glas)
ref_glas_asp=${dem_4control%.*}_${glas_file%.*}_ref_asp.csv
echo; echo $ref_glas_asp

if [ -e "$ref_glas_asp" ] ; then
    echo; echo "Run DEM co-registration to filtered ICESat-GLAS ..." ; echo
    cmd="pc_align_wrapper_3dsi.sh ${ref_glas_asp} $dem"
    echo $cmd ; echo  | tee -a $logfile
    eval $cmd  | tee -a $logfile
else
    echo; echo "No reference ICESat-GLAS points for co-reg. Exiting." ; echo
    echo; date ; echo
    exit 1    
fi

align_glas_dsm=${dem%.*}_align_glas.tif

if [ -e "$workdir/${dem_file%.*}_align/${dem_file%.*}-trans_reference-DEM.tif" ] ; then

    echo; echo "DEM co-reg to ICESat-GLAS successful" ; echo | tee -a $logfile
    ln -sf $workdir/*align/*trans_reference-DEM.tif $align_glas_dsm

    if [ $main_dir = $test_dir ] ; then
        echo; echo "Hillshade the co-reg'd DEM" ; echo
        hs_dem.sh $test_dir/$pairname/*align_glas.tif
    fi
else
    echo; echo "Failed to co-register DEM to ICESat-GLAS" ; echo | tee -a $logfile
fi

echo; date ; echo
