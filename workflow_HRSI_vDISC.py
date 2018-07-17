#!/usr/local/other/SLES11.3/miniconda/4.0.5/gcc-5.3-sp3/bin/python
###############################################
# Import and function definitions
import os, sys, osgeo, time, glob, platform, subprocess as subp
from timeit import default_timer as timer
from time import gmtime, strftime

def find_elapsed_time(start, end): # take two timer() objects and find elapsed time between them, in minutes
    elapsed_min = (end-start)/60
    return float(elapsed_min)

def run_asp(
    pairname,
    batchID,
    ASPdir,
    preLogTextFile,
    stereoDef='/discover/nobackup/projects/boreal_nga/code/stereo.default',
    searchExtList=['.ntf','.tif','.NTF','.TIF']
    ):

    #n #8 shouldnt really need any of this. Variables needed:
    # pairname; batchID; prelog text file (eh)

    # imageDir is /discover/.../ASP/batch/pairname -- imageDir is now the outDir AND inDir
    # batchDir is /discover/.../ASP/batch
    # ddir is     /discover/.../ASP/
    batchDir = os.path.join(ASPdir, 'batch{}'.format(batchID))
    imageDir = os.path.join(batchDir, pairname)
    ddir = os.path.split(ASPdir)[0] # strip off ASP to get boreal_nga
    stereoCode = os.path.join(ddir, 'code', 'evhr', 'dg_stereo.sh') # now strip off ASP (outdir name) and get code dir
    start_main = timer()
    #T:
    print ddir
    print batchDir
    print imageDir
    print stereoCode


    # hardcode stuff for now
    nodeName = platform.node()
    test = False

    # also need to read preLogText file into list
    with open(preLogTextFile, 'r') as tf:
        preLogText = tf.read()

    # For logging on the fly
    logdir = os.path.join(ddir, 'Logs')
    os.system('mkdir -p {}'.format(logdir)) # make log dir if it doesn't exist
    start_time = strftime("%Y%m%d-%H%M%S")
    lfile = os.path.join(logdir, 'batch{}__{}__{}_{}_Log.txt'.format(batchID, pairname, start_time, nodeName)) #* 2/8: putting date/time before node so it's in chrono order

    so = se = open(lfile, 'w', 0)                       # open our log file
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # re-open stdout without buffering
    os.dup2(so.fileno(), sys.stdout.fileno())           # redirect stdout and stderr to the log file opened above
    os.dup2(se.fileno(), sys.stderr.fileno())

    # print some things to the log file
    print "--LOGFILE------------------"
    print(lfile)
    print "\n"
    print "--PYTHON FILE-----------"
    print os.path.basename(__file__)
    print "\n"

    # print input parameters to log file:
    print '-------parameters:-------'
    print 'test = {}'.format(test)
    print 'batchID = {}'.format(batchID)
    print 'imageDir = {}'.format(imageDir)
    print '\nBEGIN: {}\n\n'.format(start_time)

    print "########################################\nQuery Log:\n{}\n########################################\n".format(preLogText)


    # now call the script
    print "Calling {} to perform stereo...".format(stereoCode)
    print " Parameters: {} false false {}\n\n".format(pairname, batchDir)
    #command = 'bash {} {} false false {}'.format(stereoCode, pairname, batchDir) # false for ADAPT and false for MAP -- command we've been using
    """Command above is what we've been using. Parameters below are what Paul wants to run with ANDES mini for new dg_stereo (6/1/2018):
    pairname=$1
    TEST=false
    ADAPT=false
    MAP=false
    RUN_PSTEREO=false
    batch_name='batch_andesmini'
    rpcdem=''
    NODES=false
    nodeslist=''
    SGM=false
    subpixk=25
    erode_max_size=1024

    # ANDES mini p2 parameters (6/21ish)
    test_p = 'false'
    adapt_p = 'false'
    map_p = 'false'
    runStereo_p = 'false'
    batch_p = 'batch{}'.format(batchID)
    rpc_p = ''
    nodes_p = 'false'
    nodesList_p = ''
    sgm_p = 'false'
    subpix_p = '15'
    erodeSize_p = '1024'
    """
    # 6/26 3DSI08 parameters*
    #* added parameters corrKern_p and corrTime_p (14 total parameters including pairname)
    # 3DSI08 params:
##    test_p = 'false'
##    adapt_p = 'false'
##    map_p = 'false'
##    runStereo_p = 'false'
##    batch_p = 'batch{}'.format(batchID)
##    rpc_p = ''
##    nodes_p = 'false'
##    nodesList_p = ''
##    sgm_p = 'false'
##    subpix_p = '7'
##    erodeSize_p = '0'
##    corrKern_p = '21'
##    corrTime_p = '300'

    # 7/17/2018: ANDES p3 / FL batch params
    test_p = 'false'
    adapt_p = 'false'
    map_p = 'false'
    runStereo_p = 'false'
    batch_p = 'batch{}'.format(batchID)
    rpc_p = ''
    nodes_p = 'false'
    nodesList_p = ''
    sgm_p = 'false'
    subpix_p = '21'
    erodeSize_p = '1024'
    corrKern_p = '21'
    corrTime_p = '600'

    command = 'bash {} {} {} {} {} {} {} "{}" {} "{}" {} {} {} {} {}'.format(stereoCode, pairname, test_p, adapt_p, map_p, runStereo_p, batch_p, rpc_p, nodes_p, nodesList_p, sgm_p, subpix_p, erodeSize_p, corrKern_p, corrTime_p) # parameters specified by Paul specifically for andes mini batch. temporary most likely
    #subp.check_output([command])
    print 'New command: {}\n\n'.format(command)
    os.system(command) # try this for now



    # also move the slurm.out file(s) to outASP/slurmOuts/batchID and rename them to the pairname_slurm.out
    # loop through all slurm files in pairname directory and rename/copy them to outSlurm dir -- if we are rerunning a pair process it will just name/recopy the slurm.out files to outSlurm
    inSlurmGlob = glob.glob(os.path.join(imageDir, 'slurm*out')) # list of all slurm files in pairname dir
##    print inSlurmGlob #T
    outSlurmDir = os.path.join(ddir, 'outSlurm') # doing just one big dir for outSlurm now. files will have batch names in them though
    os.system('mkdir -p {}'.format(outSlurmDir))
    print ''
    for inSlurm in inSlurmGlob: # loop through slurm files
        outSlurm = os.path.join(outSlurmDir, os.path.basename(inSlurm).replace("slurm", "batch{}__{}__slurm".format(batchID, pairname)))
        print " Copying outSlurm file {} to a new name {}".format(inSlurm, outSlurm)
        os.system('cp {} {}'.format(inSlurm, outSlurm)) # review this after we are sure it works

    print "\n Adding pair {} to completedPairs text file and recording run time information to spreadsheet".format(pairname)

    #* check for the final ovr file and if it exists, add pair to list
    finalFile = os.path.join(imageDir, '{}_ortho.tif.ovr'.format(pairname))
##    print finalFile #T
    if os.path.isfile(finalFile):
        comp_pair_dir = os.path.join(ddir, 'batchSummary')
        os.system('mkdir -p {}'.format(comp_pair_dir))
        completed_pairs_txt = os.path.join(comp_pair_dir, 'batch{}_completedPairs.txt'.format(batchID))
        with open (completed_pairs_txt, 'a') as cp:
            cp.write('{}\n'.format(pairname))
    else: print " Final ovr file ({}) does not exist. Something went wrong, please check the log.".format(finalFile)
    end_main = timer()
    total_time = find_elapsed_time(start_main, end_main)

    # add some info to the run_times csv for NCCS
    # then print batchID, pairname, total_time (minutes and hours) to csv
##    strip1size = round(os.path.getsize(fullPathStrips[0])/1024.0/1024/1024, 3) #* PC_tif here instead
##    strip2size = round(os.path.getsize(fullPathStrips[1])/1024.0/1024/1024, 3)
    # get the size of the out-PC file:
    outPC = os.path.join(imageDir, 'out-PC.tif')
    PCsize_GB = round(os.path.getsize(outPC)/1024.0/1024/1024, 3)

    run_times_csv = os.path.join(ddir, 'run_times.csv')
    with open(run_times_csv, 'a') as rt:
        rt.write('{}, {}, {}, {}, {}, {}\n'.format(batchID, pairname, total_time, (total_time/60), PCsize_GB, nodeName))

    print("\n\n-----------------------------")
    print("\n\t ")
    print("Finished processing {}".format(pairname))
    print 'End time: {}'.format(strftime("%Y%m%d-%H%M%S"))
    print "Elapsed time = {} minutes".format(round(total_time, 3))
    print("\n\t ")
    print("-----------------------------\n\n")


##    # try to close the out/err files-- http://stackoverflow.com/questions/7955138/addressing-sys-excepthook-error-in-bash-script
##    try:
##        sys.stdout.close()
##    except:
##        pass
##    try:
##        sys.stderr.close()
##    except:
##        pass


if __name__ == "__main__":
    #import sys
    # get variables being passed along from query_db and run_asp with them
    #run_asp( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8]  )
    # TESTING passing vars thru
    run_asp( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])#, sys.argv[5], sys.argv[6],  sys.argv[7], sys.argv[8],  sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12], sys.argv[13] ) # 13 arguments (plus python script)



