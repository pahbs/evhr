#!/bin/bash

batch=$1 # only command line argument is batch

logfile="/att/nobackup/mwooten3/AIST/TTE/queryLogs/batch${batch}_ADAPT_query_log.txt"

printf "\nCopying archive to DISCOVER with rsync. Output can be found in:\n$logfile\n\n"

exec >> $logfile 2>&1 # redirect standard out and error to log file

SECONDS=0


printf "\n-------------------------------------------"
printf "\nCopying archive to DISCOVER with rsync:\n\n"
printf "START: "
#start_date=date
date
#printf "\n"
cmd="nohup rsync -avxH --progress /att/nobackup/mwooten3/AIST/TTE/ASP/batch$batch-archive.tar.gz discover.nccs.nasa.gov:/discover/nobackup/projects/boreal_nga/ASP"
printf " $cmd\n\n"
eval $cmd
printf "\n\nEND: "
#end_date=date
#printf "$end_date"
date
duration=$SECONDS
printf "\nElapsed time to move batch $batch archive = $(($duration/60)) minutes\n\n"
