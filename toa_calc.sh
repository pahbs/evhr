#! /bin/bash

#Convert image DN to top-of-atmosphere reflectance
# modified from toa.sh

#Input directory containing DG Ids
dir=$1
res=4
cd $dir

#Use the image closest to nadir
#This should work for directories containing one or two IDs
nadir_id=$(nadir_id.sh . | awk '{print $1}')
#For an individual image
#nadir_id=$1
echo "Nadir: $nadir_id"

# edit here to get med res ortho
#PMedit: ortho=$(ls *${nadir_id}*_ortho_*m.tif | sort -n | head -1)
ortho=$(ls *${nadir_id}*_ortho.tif | sort -n | tail -1)
img=$ortho
echo "Image: $img"

if [ ! -e ${img%.*}_${res}m.tif ] ; then
    echo "Generating lowres image: ${img%.*}_${res}m.tif"
    gdalwarp -overwrite -r average -dstnodata 0 -tr $res $res $ortho ${ortho%.*}_${res}m.tif 
fi
img=${img%.*}_${res}m.tif

if [ -e ${img%.*}_toa.tif ] ; then
    echo "Found existing toa image: ${img%.*}_toa.tif"
    exit
fi

#Create list of available xml
#Note: older versions of dg_mosaic didn't output average MEANSUNEL, limit to subscene xml
#xml=$(ls *${nadir_id}*.xml | grep P1BS | head -1)
#dg_mosaic now writes out ABSCALFACTOR and EFFECTIVEBANDWIDTH to r100.xml
xml_list=$(ls *${nadir_id}*.xml)

#Check that we have necessary constants from xml
outxml=''
for xml in $xml_list
do
    if grep -q ABSCALFACTOR $xml && grep -q EFFECTIVEBANDWIDTH $xml && grep -q MEANSUNEL $xml ; then
        outxml=$xml
        break
    fi
done
if [ -z $outxml ] ; then
    echo "Unable to find xml with ABSCALFACTOR, EFFECTIVEBANDWIDTH, and MEANSUNEL defined"
    exit
fi
echo "XML: $outxml"

band='P'
#Calculate the scaling factor
toa.py $xml $band
c=$(toa.py $xml $band | tail -n 1)

echo "Generating new TOA image using scaling factor $c"
#For 30 m image, don't need 20 threads
image_calc --output-nodata-value 0 -d float32 -c "$c*var_0" $img -o ${img%.*}_toa.tif

#If doing fullres, could scale by 1000, then output uint16
#minval=1
#maxval=2048
#maxval=1000
#TOA reflectance is scaled range is from 0-1, want to scale back to 0-2048
#image_calc -d uint16 -c "(${c}*var_0)*(${maxval}-${minval})" $img -o ${img%.*}_toa.tif
