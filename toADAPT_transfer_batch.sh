#!/bin/bash

batch=$1 # only command line argument is batch

logfile="/discover/nobackup/projects/boreal_nga/batchSummary/batch${batch}_transferToADAPT_log.txt"
batchdir="/discover/nobackup/projects/boreal_nga/ASP/batch${batch}"
ADAPTdir="/att/nobackup/mwooten3/AIST/TTE/test_rsync"
user="mwooten3"

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

#cmd="rsync -avxHR --progress --exclude '*/*r100.xml' */*txt */*.xml */*ortho*.tif */*.ovr */*out-DEM* */out-PC.tif dsclogin.sci.gsfc.nasa.gov:${ADAPTdir}"

cmd_list=''

for pairname in `ls -d [WGQ]*` ; do
   cmd=''
   cmd="rsync -avxHR --progress "
   cmd+="$pairname/*txt "
   cmd+="$pairname/*xml "
#   cmd+="$pairname/*ortho*.tif "
#   cmd+="$pairname/*.ovr "
#   cmd+="$pairname/*out-DEM* "
#   cmd+="$pairname/out-PC.tif "
   cmd+="$user@ngalogin.nccs.nasa.gov:$ADAPTdir ; "
   cmd_list+=\ \'$cmd\'
done

echo $cmd_list
eval parallel --delay 1 -verbose -j 10 ::: $cmd_list


printf " $cmd\n\n"
eval $cmd
printf "\n\nEND: "
#end_date=date
#printf "$end_date"
date
duration=$SECONDS
printf "\nElapsed time to move batch $batch archive = $(($duration/60)) minutes\n\n"
