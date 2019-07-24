#!/bin/bash
#
# Co-register a single HRSI DSM to LiDAR dataset
# 
# Takes in 4 inputs: Pairname, Output Directory, Lidar Filename, Lidar threshold value (masks values over threshold)
#
# 3 steps
# - Masking input DEM to get control surfaces (dem_control.py) Using Lidar Mask
# - PC aligning input DEM with Reference DTM
#

#source ~/anaconda3/bin/activate py2 #For Paul
source ~/miniconda2/bin/activate py2 #My code

hostN=`/bin/hostname -s`
script_name=$(basename ${0})
pairname=${1}
main_dir=${2:-''}
ref_dsm=${3:-''}
lidar_fn=${4:-''}
thresh=${5:-''}
dem_file=${6:-'out-DEM_1m.tif'}

if [ -z $main_dir ] ; then
    main_dir=$test_dir
fi


workdir=$main_dir/$pairname
mkdir -p $workdir
dem=$workdir/${dem_file}

#Pubrepo files needed: out-DEM_*.tif, *.xml; Need to needed files to NOBACKUP
ln -sf /att/pubrepo/DEM/hrsi_dsm/v2/$pairname/out-DEM*m.tif $main_dir/$pairname
ln -sf /att/pubrepo/DEM/hrsi_dsm/v2/$pairname/*ortho*tif $main_dir/$pairname
xml_fn_list=$(ls /att/pubrepo/DEM/hrsi_dsm/v2/$pairname/*.xml)
ln -sf $xml_fn_list $main_dir/$pairname

#Create Log file info
mkdir -p $main_dir/logs_coreg2lidar
logfile=$main_dir/logs_coreg2lidar/${script_name%.*}_${hostN}_${pairname}.log

echo; echo "Script call:" | tee $logfile
echo "${script_name} ${1} ${2} ${3}"| tee -a $logfile
echo "Main dir: $main_dir" | tee -a $logfile

dem_4control=${dem%1*}4m.tif
echo $dem_4control

ref_dem='/att/gpfsfs/briskfs01/ppl/pmontesa/userfs02/data/tandemx/TDM90/mos/TDM1_90m_circ_DEM.vrt'
elev_range="-10 10"


#Run DEM control to get mask from Lidar
#Has default dem_control masking included (TOA, roughness, dz, etc)
echo; echo "Get static control DEM: masks input using LiDAR dataset masked out of the 4m version of each DEM" ; echo | tee -a $logfile
cmd="dem_control_working.py $dem_4control --dilate_con 3 -filt_param $ref_dem $elev_range -lidar_fn $lidar_fn -max_thresh $thresh"
echo $cmd ; echo  | tee -a $logfile
eval $cmd  | tee -a $logfile

dem_control=${dem_4control%.*}_control.tif
ref_coreg=${ref_dsm%.*}_control.tif
mask_file=${dem_4control%.*}_finalmask.tif

#Apply mask to the reference dataset (Lidar DTM)
echo; echo "Applying Mask from dem_control.py to the reference DSM";echo | tee -a $logfile
cmd="apply_mask_lidar.py -extent mask -out_fn $ref_coreg $ref_dsm $dem_control"
echo $cmd ; echo | tee -a $logfile
eval $cmd | tee -a $logfile

#Run Coreg on the DSM with the masked Lidar Reference#
if [ -e "$ref_coreg" ] ; then
    echo; echo "Run DEM co-registration to Lidar Filtered Reference DEM ..." ; echo
    cmd="pc_align_wrapper_3dsi.sh $ref_coreg $dem"
    echo $cmd ; echo  | tee -a $logfile
    eval $cmd  | tee -a $logfile
fi

align_lidar_dsm=${dem%.*}_align_lidar.tif

#Synlink Coreg'd DSM to main dir#
if [ -e "$workdir/${dem_file%.*}_grid_align/${dem_file%.*}-trans_reference-DEM.tif" ] ; then

    echo; echo "DEM co-reg to LiDAR successful" ; echo | tee -a $logfile
    ln -sf $workdir/*align/*trans_reference-DEM.tif $align_lidar_dsm
	
	echo; echo "Hillshade the co-reg'd DEM" ; echo
	hs_dem.sh $main_dir/$pairname/*align_lidar.tif

else
    echo; echo "Failed to co-register DEM to LiDAR" ; echo | tee -a $logfile
fi

echo; date ; echo
