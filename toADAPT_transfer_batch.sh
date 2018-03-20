#!/bin/bash

batch=$1 # only command line argument is batch

logfile="/discover/nobackup/projects/boreal_nga/batchSummary/batch${batch}_transferToADAPT_log.txt"
batchdir="/discover/nobackup/projects/boreal_nga/ASP/batch${batch}"

printf "\nCopying results from DISCOVER to ADAPT with rsync. Output can be found in:\n$logfile\n\n"

exec >> $logfile 2>&1 # redirect standard out and error to log file

SECONDS=0


printf "\n-------------------------------------------"
printf "\nCopying archive to ADAPT with rsync:\n\n"
printf "START: "
#start_date=date
date
cd ${batchdir}
printf "\n Current directory:"
pwd
printf "\n"

cmd="rsync -avxHR --progress --exclude '*/*r100.xml' */*txt */*.xml */*ortho*.tif */*.ovr */*out-DEM* */out-PC.tif dsclogin.sci.gsfc.nasa.gov:/att/nobackup/mwooten3/AIST/TTE/test_rsync" #/att/pubrepo/DEM/hrsi_dsm/test_rsync"

#cd ${batchdir}
#pwd

#for d in * ; do
#    echo "$d"
#    cmd="rsync -avxH -R -p --progress ${d}/out-PC.tif dsclogin.sci.gsfc.nasa.gov:/att/nobackup/mwooten3/TTE/test_rsync"
#    echo $cmd
#    eval $cmd
#done


printf " $cmd\n\n"
eval $cmd
printf "\n\nEND: "
#end_date=date
#printf "$end_date"
date
duration=$SECONDS
printf "\nElapsed time to move batch $batch archive = $(($duration/60)) minutes\n\n"
