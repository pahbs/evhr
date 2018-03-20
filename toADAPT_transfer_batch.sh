#!/bin/bash

batch=$1 # only command line argument is batch

logfile="/discover/nobackup/projects/boreal_nga/batchSummary/batch${batch}_transferToADAPT_log.txt"
batchdir = "/discover/nobackup/projects/boreal_nga/ASP/batch${batch}"

printf "\nCopying results from DISCOVER to ADAPT with rsync. Output can be found in:\n$logfile\n\n"

exec >> $logfile 2>&1 # redirect standard out and error to log file

SECONDS=0


printf "\n-------------------------------------------"
printf "\nCopying archive to ADAPT with rsync:\n\n"
printf "START: "
#start_date=date
date
#printf "\n"
cmd="nohup rsync -avxH -R --progress ${batchdir}/./*/*txt ${batchdir}/./*/*.xml ${batchdir}/./*/*ortho*.tif ${batchdir}/./*/*.ovr ${batchdir}/./*/*out-DEM* ${batchdir}/./*/out-PC.tif dsclogin.sci.gsfc.nasa.gov:/att/nobackup/mwooten3/TTE/test_rsync"#/att/pubrepo/DEM/hrsi_dsm/test_rsync"
printf " $cmd\n\n"
eval $cmd
printf "\n\nEND: "
#end_date=date
#printf "$end_date"
date
duration=$SECONDS
printf "\nElapsed time to move batch $batch archive = $(($duration/60)) minutes\n\n"
