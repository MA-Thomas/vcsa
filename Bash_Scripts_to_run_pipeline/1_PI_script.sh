#!/bin/bash
# This is is part of the collection (1_PI_script.sh, 1_reSubmission_shell_script.sh, 1_dessaJobSubmissions_Shell_Script.sh )

# Currently, 1_reSubmission_shell_script.sh is not used and I manually restart when needed by resubmitting 1_PI_script.sh

##########################################################################################################################################
######################### This section is specific to Marcus Thomas. Please update for your compute cluster. #############################
cd /home/marcust/param_infer_122114/    # This is where all of the matlab files (including global solvers) and bash files are stored for Marcus Thomas.

 
# Define String Path Variables. If altered here, also alter in "1_dessaJobSubmissions_Shell_Script.sh".

folderName="CCMV_syn/" # CCMV virus, synthetic data
homePath="/home/marcust/"      # my home folder on the cluseter
scratchPath="/scratch/marcust/" # all users on department cluster have access to a "scratch" folder for storing temporary data. 
                                # If this is not applicable to you, use an appropriate temp folder.

scratchWorkFolder=$scratchPath$folderName  # e.g. "/scratch/marcust/CCMV_syn/"
homeWorkFolder=$homePath$folderName   # e.g. "/home/marcust/CCMV_syn/"

# The simulator (Dessa jar file "DESSA_2017_SAXS_ReducedTimes_noLikelihoodPrintouts.jar"), the template xml file 
#  ("CCMV_156_00_5capsids.xml" which is copied and altered during optimization), the CRYSOL based
# form factor for a single CCMV subunit (.mat file), and the ground truth SAXS experiment (.mat file) should be stored in:
# /home/marcust/simulator/

# If you wish to update Dessa and/or the template xml file (e.g. to simulate a system larger than 5 capsids), you must add the 
# updated file names to the appropiate structs in begin.m

# begin.m should be examined carefully, with any file or folder names specific to Marcus Thomas updated to reflect your current file system.

# Lastly, this pipeline expects there to be 3 folders in your homeWorkFolder at the start of a search.
# 1. storage/  (this should be an empty folder)
# 2. dessaFolderOnHome/  (this should be an empty folder)
# 3. err/  (this should contain 4 empty files: outFileJobs.txt, outFile.txt, errFileJobs.txt, errFile.txt)

# How to run this script:
# 1. Submit this script to cluster as a job, along with paths to output file and err file locations.

# e.g.,  sbatch -p rs1 -N 1 -n 8 --mem=12Gb -o /home/marcust/CCMV_syn_SNOBFIT/err/outFile.txt -e /home/marcust/CCMV_syn_SNOBFIT/err/errFile.txt -t 1-71:59:00 1_PI_script.sh


# 2. After a few seconds, a few files will be generated in your homeWorkFolder by begin.m/preLoop.m/inLoop1.m (called below) including the initial
#    set of random points to run simulations at, xFile.mat.
# 3. Check err/outFile.txt. When it displays "PI_script.sh: Now at the stage where dessaJobSubmissions_Shell_Script.sh takes over. ",
#    you can submit the job script "1_dessaJobSubmissions_Shell_Script.sh". This script will run the newly generated simulation jobs.
#    To alter the ground truth parameter or the number of repeated experiemts per point, see optInfo.BS or simInfo.repeat in the file begin.m
#    To alter the initially selected random points, see x in the file inLoop1.m

# 4. 1_PI_script.sh is designed to wait (see while loop below) until all the jobs have finished running before moving on to the Matlab file 
#    compObj.m which computes the RMSDs of the completed simulations. It continually checks the automatically created file jobs_running.txt (defined below). 
#    When all jobs have finished running, a variable in that file will take on the value "2". This only happens when it is detected that 
#    there is only 1 job still running on the cluster, namely, 1_PI_script.sh. 
#    THIS MAY NEED TO BE ALTERED FOR OTHER USERS.

# 5. 1_PI_script.sh will call a Matlab file named "inLoop2.m" or something similar. This file performs the optimization step and returns
#    the next points to simulate as the file xFile.mat (i.e. the next set of DESSA job submission files to be created by inLoop1.m during 
#    the next round)

# 6. If the walltime (user specified limit to how long 1_PI_script.sh can run - see maxTime below) is near, 1_PI_script.sh will not move on
#    to calling the next Matlab file in the loop. It will instead prepare to stop running. The automatically created file "resub.txt" will be 
#    updated, allowing 1_PI_script.sh to be restarted at the appropriate stage. When restarting, be sure to set the variable "wallTimeSoon" 
#    back to 0 within "reSub.txt".

# NOTE: the variable "count" specifies which round of search we are at.
#        "count" is incremented each time the loop calls the file inLoop1.m. Its value is saved in your homeWorkFolder, in 
#       the automatically created files: loop_eval_criteria.mat and loop_eval_criteria.txt. 
#       If you need to restart 1_PI_script.sh in inLoop1.m, be sure that "count" is stored as the 1 less than the current round.


# NOTE: To restart the search in any previous round, see the automatically created folder "storage/" for the needed files from the relevant
#       round. These files should be (after appropriate renaming) put back in homeWorkFolder.

#########################################################################################################################################
#########################################################################################################################################


# Marcus Thomas
# Define Shell Variables and Functions.
startTime=$(date +%s)           # num seconds since 1970-01-01 00:00:00 UTC
maxTime=`expr 86400 - 60`  # corresponds to 23:59:00 below (60s less than 24hours)
currTime() {                             # function to be queried throughout this script.
  date +"%s"
}                                       
minTimeLeft=240                              #when there's only $minTimeLeft seconds until walltime is reached, prepare to resubmit Parameter Inference 
timeNeededForCompObj=7200      #at most, this much time needed for compObj.m
timeNeededForInLoop2=3600       #at most, this much time needed for InLoop2.m (to be used with inLoop1.m check at end of loop too)
sleepTime=30                                    #time to wait before checking if jobs finished 1) being submitted and, 2) running
queryQueue=squeue                         #qstat -Q in PBS
clusterStatus=sinfo                           #qstat -a in PBS

walltimeSoon=0             #0 indicates no, 1 indicates yes. This variable (in $reSubStatusFileOnHome) will be queried by
                                         #[reSubmission_Shell_Script.txt] which lives on the head node.
allDone=0                        #0 indicates no, 1 indicates yes. When EVERYTHING is done, set this to 1 to alert [reSubmission_Shell_Script.txt].


echo "DON'T RUN OTHER JOBS ON THIS QUEUE AS THAT WILL RENDER QUERIES UNMEANINGFUL"
echo
echo "#### in this script signifies code related to handling case where walltime limit is reached"
echo
echo "STARTING"
echo



fileInfoPrefix=${folderName%/}  #Remove the suffix starting with "/" using the % operator, e.g. "CCMV_syn"


dessaFolder="dessaFolderOnHome/"
#errFolder="err/"                                                   # Must already exist on home workfolder, and contain errFile.txt and outFile.txt
loopFile="loop_eval_criteria.txt"

jlist="job_list.txt"
jobListFileOnHome=$homePath$folderName$jlist   #"/home/marcust/CCMV_syn/job_list.txt"
##jobListFileOnHome="/home/marcust/HPV_snobfit_SAXS/job_list.txt"

jobsRunningFileOnHome="jobs_running.txt"    #running=(0,1,2) 
                                                                                 #0-none submitted, 1-all submitted, 2-all finished running
star="*"    #In case I want to operate on contents of a dir, but not the dir itself

# Make sure error folder exists
#[ -d $homeWorkFolder$errFolder ] || mkdir $homeWorkFolder$errFolder

#### File Paths For Case Where Script Not Completed by the walltime limit.
reSubStatusFileOnHome="resub.txt"         # This file, which lives on /home/, contains either reSub=0 (first run case) or reSub=1 (post first run)
                                                                        # If reSub=1, we must complete unfinished parts of this script
                                                                        # Also contains reStartInLoop1/reStartCompObj/reStartInLoop2 which are ONLY to be considered 
                                                                        # if reSub=1.
                                                                        # If reSub=1 AND reStartInLoop1=1, then let inLoop1.m me called. 
                                                                        # Same for reStartCompObj and reStartInLoop2.m and cleanup.m.

#### If [resub.txt] already exists do nothing, otherwise this is the first run case, so create it.
if [ -f $homeWorkFolder$reSubStatusFileOnHome ]; then echo;
else
    echo reSub=0 > $homeWorkFolder$reSubStatusFileOnHome
    echo reStartInLoop1=0 >> $homeWorkFolder$reSubStatusFileOnHome
    echo reStartCompObj=0 >> $homeWorkFolder$reSubStatusFileOnHome
    echo reStartInLoop2=0 >> $homeWorkFolder$reSubStatusFileOnHome
    echo walltimeSoon=0 >> $homeWorkFolder$reSubStatusFileOnHome
    echo cleanup=0 >> $homeWorkFolder$reSubStatusFileOnHome
    echo allDone=0 >> $homeWorkFolder$reSubStatusFileOnHome
fi

#### Set shell loop variables to values in [resub.txt]
source $homeWorkFolder$reSubStatusFileOnHome
reSubmitted=$reSub
runInLoop1=$reStartInLoop1
runCompObj=$reStartCompObj
runInLoop2=$reStartInLoop2
runCleanup=$cleanup

# Either [jobs_running.txt] already exists, OR this is the first run case, so create it.
[ -f $homeWorkFolder$jobsRunningFileOnHome ] || echo running=0 > $homeWorkFolder$jobsRunningFileOnHome

# Make sure the workfolder path assigned in begin.m is identical to this one.
# Only if it exists, delete old /scratch/ workfolder.
[ -d  $scratchWorkFolder ] && rm -r $scratchWorkFolder

# Either /scratch/marcust/ already exists OR create it.
[ -d /scratch/marcust/ ] || mkdir /scratch/marcust/

# Now make the subfolder.
mkdir $scratchWorkFolder

#### TODO: not urgent. give all file paths to begin, store in fileInfo struct so I don't need to pass arguments to matlab functions here.
if [ $reSubmitted -eq 0 ]; then
    /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try begin('$homeWorkFolder','$scratchWorkFolder','$fileInfoPrefix'); catch; end; quit"
    /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try preLoop('$homeWorkFolder'); catch; end; quit"

    ## May12,2016 update: preLoop.m, inLoop1.m now take homeWorkFolder, not scratchWorkFolder.
    ## This made sense because inLoop1.m kept updating minScore from /scratch/ (minScore=inf),
    ## while it should have updated after initial round from /home/. Additionally, because there 
    ## isn't much saving going on it makes sense to just utilize /home/ and leave /scratch/ for
    ## prepareJob.m which can access it through fileInfo.workFolder.
    # If begin.m and preLoop.m are run:
    # (Update: Don't delete from /scratch/ yet since inLoop1.m needs stuff on /scratch/) 
    # This will result in $homeWorkFolder/ as the path prefix
    # From begin.m and preLoop.m, this includes paramFile.mat, xFile.mat, loop_eval_criteria.txt, paramFile.mat

    #[ -d  $scratchWorkFolder ] && [ "$(ls $scratchWorkFolder )" ] && cp -r $scratchWorkFolder$star $homeWorkFolder


    echo "PI_script.sh: begin.m and preLoop.m finished"
    noImprov=0;
    optInfo_maxIter=10;
    count_iterations=0;

elif [ $reSubmitted -eq 1 ]; then
    source $homeWorkFolder$loopFile
    noImprov=$noImprove
    optInfo_maxIter=10;
    count_iterations=$count
else
    echo
    #ERROR HANDLING GOES HERE. $reSubmitted should ALWAYS be 0 or 1.
fi



# echo "Jan 22 evening - Test Complete!"
# exit
# : <<'COMMENTED'
 optInfo_maxIter=30;
# BEGIN MAIN LOOP  
###############################################################
###############################################################
#while [ $noImprov -le $optInfo_maxIter ]; do         #For use with inLoop2.m
while [ $count_iterations -le $optInfo_maxIter ]; do #For use with inLoop2_BayesianOpt.m
#while [ $noImprov -lt 3 ]; do

        # If inLoop1.m is run:
        # (Matlab system command: Copies all contents of scratchWorkFolder to /home/, then here,  delete contents of scratchWorkFolder) 
        # This will result in $homeWorkFolder/ as the path prefix
        # From inLoop1.m, this includes job_list.txt, dataList.mat, jobNum.mat, iterPrefix.mat, loop_eval_criteria.txt

        if [ $reSubmitted -eq 1 ]; then
            if [ $runInLoop1 -eq 1 ]; then                
                /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try inLoop1('$homeWorkFolder'); catch; end; quit"
                #[ -d  $scratchWorkFolder ] && [ "$(ls $scratchWorkFolder )" ] && cp -r $scratchWorkFolder$star $homeWorkFolder
                
                echo "PI_script.sh: inLoop1.m Finished (just resubmitted - PI_script.sh)"
                echo
                sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                echo reSub=0 >> $homeWorkFolder$reSubStatusFileOnHome                    
                sed -i '/reStartInLoop1/d' $homeWorkFolder$reSubStatusFileOnHome     
                echo reStartInLoop1=0 >> $homeWorkFolder$reSubStatusFileOnHome
                source $homeWorkFolder$reSubStatusFileOnHome
                reSubmitted=$reSub
                runInLoop1=$reStartInLoop1

            fi
        else      
            /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try inLoop1('$homeWorkFolder'); catch; end; quit"
            #[ -d  $scratchWorkFolder ] && [ "$(ls $scratchWorkFolder )" ] && cp -r $scratchWorkFolder$star $homeWorkFolder
            
            echo "PI_script.sh: inLoop1.m Finished (non-resubmitted - PI_script.sh)"
            echo
        fi

        
        # *************************************************************************************
        # ************  [dessaJobSubmissions_Shell_Script.sh]  SUBMITS DESSA JOBS here ***************
        # /scratch/ not used further in this script, aside from the individual dessa job submissions.
        echo "PI_script.sh: Now at the stage where dessaJobSubmissions_Shell_Script.sh takes over. "
        echo       

        # Set value of 'running'. 
        # Remember, 0=none or some submitted, 1=all submitted, 2=all finished running. Set in [dessaJobSubmissions_Shell_Script.sh]
        source $homeWorkFolder$jobsRunningFileOnHome
        while [ $running -eq 0 ]; do            
            sleep $sleepTime
            sleep $sleepTime
            sleep $sleepTime
            sleep $sleepTime

            #### Perform time checks while waiting for jobs to finish running. 
            cTime="$(currTime)"
            elapsedTime=`expr $cTime - $startTime`
            timeLeft=`expr $maxTime - $elapsedTime`
            if [ $timeLeft -le $minTimeLeft ]; then          
                sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
                echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome     # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

                sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome

                sleep $minTimeLeft
            fi
            ####

            source $homeWorkFolder$jobsRunningFileOnHome
            val="PI_script.sh: Waiting for dessaJobSubmissions_Shell_Script.sh to notify us of changes to running variable. It's current value is "
            echo $val$running
            echo            
        done      
        
        echo "PI_script.sh: running is no longer 0"
        echo
         
        source $homeWorkFolder$jobsRunningFileOnHome
        if [ $running -ne 1 ]; then
            if [ $running -ne 2 ]; then
                echo "PI_script.sh: Problem submitting or running all jobs. Exiting." >> $homeWorkFolder$jobsRunningFileOnHome
                exit
            fi
        fi
        echo "PI_script.sh: No Problem submitting jobs."
        echo  
        while [ $running -ne 2 ]; do
        	    
        	    echo "PI_script.sh: In while loop which waits for submitted jobs to finish running (running=2). In PI_script.sh."
        	    echo "PI_script.sh: Currently, running is: "
        	    echo $running
        	    echo
        	    source $homeWorkFolder$jobsRunningFileOnHome
        	    sleep 60

                cTime="$(currTime)"
                elapsedTime=`expr $cTime - $startTime`
                timeLeft=`expr $maxTime - $elapsedTime`

                #### If walltime limit is near, save current state on /home/, 
                #### resubmit PI script (walltimeSoon), let code sleep till walltime reached
                if [ $timeLeft -le $minTimeLeft ]; then          
                    sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
                    echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome              # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

                    sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                    echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome

                    sleep $minTimeLeft
                fi  
        done

        source $homeWorkFolder$jobsRunningFileOnHome
        echo "PI_script.sh: Running is now: "
        echo $running
        echo
        #### Perform another time check in case prev while loop was bypassed due to running=2.
        cTime="$(currTime)"
        elapsedTime=`expr $cTime - $startTime`
        timeLeft=`expr $maxTime - $elapsedTime`
        if [ $timeLeft -le $minTimeLeft ]; then          
            sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
            echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome              # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

            sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
            echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome

            sleep $minTimeLeft
        fi

        echo "PI_script.sh: ALL JOBS for this round have been submitted and finished running."
        echo "PI_script.sh:  [job_list.txt] has been deleted from home directory by dessaJobSubmissions_Shell_Script.sh."
        echo
                   
            


        # TODO (done!): 
        #       cleanup.m call (from within inLoop2.m and snobopt.m) commented out. It generated the tgz storage folders.

        ## UPDATE FOR THIS BLOCK Jan.26 2016: DONE BY JOB FILE ITSELF (prepareJob.m)
        # COPY DESSA DATA FROM /scratch/ TO /home/. (copy contents of scratchWorkFolder, then delete contents)
        # If reSubmission was initiated to complete PI_script.sh starting from compObj.m, these operations will have just been done.
        # (No harm in running them again)
        # This will result in the data path prefixes $homeWorkFolder$dessaFolderOnHome
        #[ -d  $scratchWorkFolder ] &&  [ "$(ls $scratchWorkFolder )" ] && cp -r $scratchWorkFolder$star $homeWorkFolder$dessaFolderOnHome
        #[ -d  $scratchWorkFolder ] &&  [ "$(ls $scratchWorkFolder )" ] && rm -r $scratchWorkFolder

        #echo "Exiting NOW. REMOVE THIS LATER"
        #exit 

        #                                                                       ********** compObj.m *********** 
        #### Perform another time check. This is where I need to make a judgement call for the max possible runtime for compObj.m
        cTime="$(currTime)"
        elapsedTime=`expr $cTime - $startTime`
        timeLeft=`expr $maxTime - $elapsedTime`
        if [ $timeLeft -le $timeNeededForCompObj ]; then          
            sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
            echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome              # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

            sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
            echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome         
                    
            sed -i '/reStartCompObj/d' $homeWorkFolder$reSubStatusFileOnHome     
            echo reStartCompObj=1 >> $homeWorkFolder$reSubStatusFileOnHome                        
            sleep $sleepTime
            sleep $sleepTime
            sleep $sleepTime
            exit
        fi
        
        # If Applicable, run compObj.m
        if [ $reSubmitted -eq 1 ]; then
            if [ $runCompObj -eq 1 ]; then
                # compObj.m  - saves a .mat file containing y
                /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try compObj('$homeWorkFolder'); catch; end; quit"                
                echo "PI_script.sh: compObj.m Finished (just resubmitted - PI_script.sh)"
                echo
                sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                echo reSub=0 >> $homeWorkFolder$reSubStatusFileOnHome                   
                sed -i '/reStartCompObj/d' $homeWorkFolder$reSubStatusFileOnHome     
                echo reStartCompObj=0 >> $homeWorkFolder$reSubStatusFileOnHome 

                source $homeWorkFolder$reSubStatusFileOnHome
                reSubmitted=$reSub
                runCompObj=$reStartCompObj    
            fi
        else      
            /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try compObj('$homeWorkFolder'); catch; end; quit"
            echo "PI_script.sh: compObj.m Finished (non-resubmitted - PI_script.sh)"
            echo
        fi


        
        #                                                                       ********** InLoop2.m ***********
        #### Perform another time check. This is where I need to make a judgement call for the max possible runtime for InLoop2.m
        cTime="$(currTime)"
        elapsedTime=`expr $cTime - $startTime`
        timeLeft=`expr $maxTime - $elapsedTime`
        if [ $timeLeft -le $timeNeededForInLoop2 ]; then          
            sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
            echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome              # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

            sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
            echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome         
                    
            sed -i '/reStartInLoop2/d' $homeWorkFolder$reSubStatusFileOnHome     
            echo reStartInLoop2=1 >> $homeWorkFolder$reSubStatusFileOnHome                        
            sleep $sleepTime
            sleep $sleepTime
            sleep $sleepTime
            exit
        fi
        
        

        # If Applicable, run inLoop2.m
        if [ $reSubmitted -eq 1 ]; then
            if [ $runInLoop2 -eq 1 ]; then
                # inLoop2.m  - calls snobfit, saves loop evaluation criteria to a .mat file [loop_eval_criteria.txt]
                /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try inLoop2('$homeWorkFolder'); catch; end; quit"                
                echo "PI_script.sh: inLoop2_BayesianOpt.m Finished (just resubmitted - PI_script.sh)"
                echo
                sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                echo reSub=0 >> $homeWorkFolder$reSubStatusFileOnHome                   
                sed -i '/reStartInLoop2/d' $homeWorkFolder$reSubStatusFileOnHome     
                echo reStartInLoop2=0 >> $homeWorkFolder$reSubStatusFileOnHome

                source $homeWorkFolder$reSubStatusFileOnHome
                reSubmitted=$reSub
                runInLoop2=$reStartInLoop2 
            fi
        else
            echo "PI_script.sh: About to run inLoop2.m. " #" or exit if y.mat does not exist."   
            #[ -f $homeWorkFolder$ymat ] || exit
            /opt/matlab/9.2/bin/matlab -nodisplay -nosplash -r "try inLoop2('$homeWorkFolder'); catch; end; quit"
            echo "PI_script.sh: inLoop2_BayesianOpt.m Finished (non-resubmitted - PI_script.sh)"
            echo
        fi

        #echo "Exiting NOW. REMOVE THIS LATER"
        #exit

        # Set shell loop variables to values in [loop_eval_criteria.txt]
        source $homeWorkFolder$loopFile
        noImprov=$noImprove

        echo "PI_script.sh: noImprove is: "
        echo $noImprov
        echo


        #                                                                       ********** End of Loop ***********
        #### Perform another time check. Might need to restart with inLoop1.m
        cTime="$(currTime)"
        elapsedTime=`expr $cTime - $startTime`
        timeLeft=`expr $maxTime - $elapsedTime`
        if [ $timeLeft -le $timeNeededForInLoop2 ]; then          
            sed -i '/walltimeSoon/d' $homeWorkFolder$reSubStatusFileOnHome     # delete walltimeSoon variable from file
            echo walltimeSoon=1 >> $homeWorkFolder$reSubStatusFileOnHome       # alerts [reSubmission_shell_script.txt] to submit PI_script.sh

            sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
            echo reSub=1 >> $homeWorkFolder$reSubStatusFileOnHome         
                    
            sed -i '/reStartInLoop1/d' $homeWorkFolder$reSubStatusFileOnHome     
            echo reStartInLoop1=1 >> $homeWorkFolder$reSubStatusFileOnHome                        
            sleep $sleepTime
            sleep $sleepTime
            sleep $sleepTime
            exit
        fi
        
        # Cleanup. Put dessa .dat files in archive (.tgz) folder stored in homeWorkFolder. Empty dessaFolder for use in next round.
        # Sometimes we may wish to start this script from 'cleanup'
        if [ $reSubmitted -eq 1 ]; then         
            if [ $runCleanup -eq 1 ]; then
                fname="data_round_"
                fextension=".tar.gz"
                concat_name="$fname$count$fextension"
                archiveName=$homeWorkFolder$concat_name
                tar czf $archiveName $homeWorkFolder$dessaFolder
                #rm -r $homeWorkFolder$dessaFolder$star

 
                ### Reset running variable to 0 for the next iteration.
                sed -i '/running/d' $homeWorkFolder$jobsRunningFileOnHome
                echo running=0 >> $homeWorkFolder$jobsRunningFileOnHome

                sed -i '/reSub/d' $homeWorkFolder$reSubStatusFileOnHome                   
                echo reSub=0 >> $homeWorkFolder$reSubStatusFileOnHome                    
                sed -i '/cleanup/d' $homeWorkFolder$reSubStatusFileOnHome     
                echo cleanup=0 >> $homeWorkFolder$reSubStatusFileOnHome 
                #noImprov=3
                #echo "PI_script.sh: noImprov is 3, so we're done! "
            fi
        else
            fname="data_round_"
            fextension=".tar.gz"
            concat_name="$fname$count$fextension"
            archiveName=$homeWorkFolder$concat_name
            tar czf $archiveName $homeWorkFolder$dessaFolder
            #rm -r $homeWorkFolder$dessaFolder$star

 
            ### Reset running variable to 0 for the next iteration.
            sed -i '/running/d' $homeWorkFolder$jobsRunningFileOnHome
            echo running=0 >> $homeWorkFolder$jobsRunningFileOnHome
            #noImprov=3
            #echo "PI_script.sh: noImprov is 3, so we're done! "
        fi     
               
        done
        # END MAIN LOOP 
        ###############################################################
        ###############################################################

sed -i '/allDone/d' $homeWorkFolder$reSubStatusFileOnHome     # delete allDone variable from file
echo allDone=1 >> $homeWorkFolder$reSubStatusFileOnHome              # alerts [reSubmission_shell_script.txt] that everthing is done

exit    

#COMMENTED


# ----------------------------------------------------------------------------------------------------------
# GENERAL COMMENTS ABOUT RUNNING THIS SCRIPT


# begin_preLoop.sh will run very quickly (save paramFile, generate parameter grid)
# inLoop1.sh will also run quickly       (saves job_list.txt and data_list.txt)


