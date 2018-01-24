#!/bin/bash

# Part of the collection containing (1_PI_script.sh, 1_reSubmission_shell_script.txt, 1_dessaJobSubmissions_Shell_Script.txt)

# This script is run from the head node. 
# It queries the existence of the file $jobListFileOnHome created in 1_PI_script.sh.


: <<'COMMENTED'
echo "dessaJobSubmissions... running "
exit
COMMENTED

homeWorkFolder="/home/marcust/CCMV_syn/"


jobListFileOnHome="job_list.txt"
jobsRunningFileOnHome="jobs_running.txt"
reSubStatusFileOnHome="resub.txt"
longTermStorage="longTermStorage/"
errFileJobs="err/errFileJobs.txt"                                    # This file must already exist (created manually in homeWorkFolder)
outFileJobs="err/outFileJobs.txt"                                    # This file must already exist (created manually in homeWorkFolder)
dessaFolder="dessaFolderOnHome/"

sleepTimeJob=20               # Time to wait before checking [job_list.txt] existence, also time to wait before checking if all jobs completed.
sleepTimeQueue=20          # Time to wait for num of jobs hanging to go below $maxOnQueue
maxOnQueue=120               #max number of jobs hanging on queue
socketPause=30                 #time to wait after each job submission so as not to have too many sockets open. NOT NEEDED with SLURM.
jobfileCounter=0                # Number of times it finds&submits the jobs_list.txt file


# Either long term storage folder exists on /home/, OR create it
[ -d $homeWorkFolder$longTermStorage ] || mkdir $homeWorkFolder$longTermStorage 



allDone=0
while [ $allDone -eq 0 ]; do
  
  # Set allDone
  if [ -f $homeWorkFolder$reSubStatusFileOnHome ]; then
  echo
  echo
  source $homeWorkFolder$reSubStatusFileOnHome
    
  if [ $allDone -eq 0 ]; then
    
    if [ -f $homeWorkFolder$jobListFileOnHome ]; then
          echo "job_list.txt exists!"
          echo 
          ((jobfileCounter++))
          # Save $homeWorkFolder$jobListFileOnHome to a read only folder for long term storage since it will be deleted by PI_script.sh
          cp $homeWorkFolder$jobListFileOnHome $homeWorkFolder$longTermStorage 


          # SUBMIT JOBS [filenames in jobList file] *****************************************
          #                                                                     *****************************************
          
          linesRead=0                                                                             #total number of jobs already submitted to queue from [job_list.txt]  
          unset prevParamPoint
          unset currParamPoint
          while read -r line; do
              echo "In while loop which iterates over files in job_list.txt"
              echo                

              # see: http://stackoverflow.com/questions/19482123/extract-part-of-a-string-using-bash-cut-split
              line2=${line##$homeWorkFolder}                                               # remove string prefix "/home/..."
              line3=${line2//"/"/_}                                                                     # replace all instances of "/" with "_" in variable line2
              line4=${line3/".sh"/".dat"}                                                           # replace single instance of ".sh" with ".dat" in line3
              echo $homeWorkFolder$dessaFolder$line4
              echo $line
              echo
              [ -f $homeWorkFolder$dessaFolder$line4 ] && continue       # if job was already run and .dat file is seen in dessaFolderOnHome, continue.

              fname="data_round_"
              round=$(echo $line4 | grep -o -P '(?<=iter).*(?=_grid)') # The round is the same whether or not the point changes.
              pname="_point_"
              Point=$(echo $line4 | grep -o -P '(?<=grid).*(?=_)')
              Point=$(echo $Point | grep -o -P '(?<=).*(?=_)') #See variable, prevParamPoint, below for how this code works.
              fextension=".tar.gz"
              concat_name="$fname$round$pname$Point$fextension"
              [ -f $homeWorkFolder$dessaFolder$concat_name ] && continue   # if job is part of previously completed point (.tar.gz of all repeats exists), continue.     


              counter=$(squeue | grep -c marcust)                                         #number of jobs currently in the queue
              while [ $counter -ge $maxOnQueue ]; do
                  sleep $sleepTimeQueue 
                  counter=$(squeue | grep -c marcust)
              done            
              
              # ------------------------------------------------------------------------------------------------------------------------
              # This is where we compress all .dat files (repeated simulations) at the previous parameter point - if appropriate.
              # If currParamPoint != prevParamPoint, we should pause here (to give all simulations at prevParamPoint time to complete),
              # then compress the corresponding .dat files into a single .tgz file in dessaFolder.
              x="temp_string" # See https://stackoverflow.com/questions/3601515/how-to-check-if-a-variable-is-set-in-bash
              if [ -z ${prevParamPoint+x} ]; then 
                echo "prevParamPoint is unset"; 
                
                # The parameter point and round are given in the simulation job filename, i.e. variable 'line4'
                # Set prevParamPoint as substring between "grid" and the last "_"
                prevParamPoint=$(echo $line4 | grep -o -P '(?<=grid).*(?=_)')
                # Only keep the part of the substring before "_"
                prevParamPoint=$(echo $prevParamPoint | grep -o -P '(?<=).*(?=_)')
                repeatNum=1 #First repeat at first param point about to be submitted.
                echo $repeatNum
              else 
                echo "prevParamPoint is set to '$prevParamPoint'"

                # Eval currParamPoint
                currParamPoint=$(echo $line4 | grep -o -P '(?<=grid).*(?=_)')
                # Only keep the part of the substring before "_"
                currParamPoint=$(echo $currParamPoint | grep -o -P '(?<=).*(?=_)')
                echo $currParamPoint

                # If currParamPoint != prevParamPoint, enter while loop to bide time, once all .dat files are there, compress them.
                #       Then exit while loop. 
                if [ "$currParamPoint" -ne "$prevParamPoint" ]; then
                  echo "currParamPoint and prevParamPoint are different. Simulations at new point about to start after outputs at prev point are compressed to a .tgz"
                  stopWaiting=0
                  while [ $stopWaiting -eq 0 ]; do
                    sleep 90
                    
                    #Value of repeatNum is now equal to the max number of repeated simulations per point.
                    # Check to see if the number of .dat files in dessaFolder equals repeatNum.
                    numDatFiles=$(ls -l $homeWorkFolder$dessaFolder*.dat | wc -l) 
                    echo "numDatFiles is:"
                    echo $numDatFiles
                    if [ "$numDatFiles" -ge "$repeatNum" ]; then                      
                      fname="data_round_"
                      round=$(echo $line4 | grep -o -P '(?<=iter).*(?=_grid)') # The round is the same whether or not the point changes.
                      pname="_point_"
                      fextension=".tar.gz"
                      concat_name="$fname$round$pname$prevParamPoint$fextension"
                      

                      currDir=$(pwd)
                      cd $homeWorkFolder$dessaFolder 
                      tar czf $homeWorkFolder$dessaFolder$concat_name *.dat
                      rm *.dat # This should remove .dat files only from within this folder.
                      cd $currDir 
                      # TO Extract (within say, compObj.m), the code the same, except 'c' is replaced with 'x' for extract.
                      # currDir=$(pwd)
                      # cd $homeWorkFolder$dessaFolder 
                      # tar xzf $homeWorkFolder$dessaFolder$concat_name *.dat
                      # cd $currDir

                      stopWaiting=1;
                    fi                    
                    
                  done

                  # Update prevParamPoint and repeatNum
                  prevParamPoint=$currParamPoint
                  repeatNum=1 #First repeat at current param point about to be submitted. 

                else
                  
                  # Update prevParamPoint and repeatNum
                  prevParamPoint=$currParamPoint
                  repeatNum=$((repeatNum+1)) #Next repeat 
                fi
                
              fi
              echo
              echo
              # ------------------------------------------------------------------------------------------------------------------------

              # Whatever memory (e.g. 8G) allocated to SLURM here must be at least 4GB larger than what's allocated
              #  to java (fileInfo.command in begin.m)
              sbatch -p rs1 --mem=12G -t 1-23:59:00 -o $homeWorkFolder$outFileJobs -e $homeWorkFolder$errFileJobs $line        
              ((counter++))
              ((linesRead++))
              echo "another job just submitted..."
              echo        
         
          done < "$homeWorkFolder$jobListFileOnHome"
          # JOBS SUBMITTED [filenames in jobList] *******************************************
          #                                                                     ******************************************* 




         # WAIT FOR JOBS TO FINISH RUNNING **********************************
         #                                                                 **********************************
          echo "Finished job submissions. Now waiting for them to finish running. "
          echo
          sed -i '/running/d' $homeWorkFolder$jobsRunningFileOnHome             # delete 'running' variable from file
          echo running=1 >> $homeWorkFolder$jobsRunningFileOnHome             # alerts [PI_script] that all jobs have been submitted
          source $homeWorkFolder$jobsRunningFileOnHome

          while [ $running -eq 1 ]; do 
              echo "In while loop that queries the variable 'running'. Its value is: "
              echo $running
              echo
              sleep $sleepTimeJob

              # If jobs still running on rs1. PI_script.sh should always be running, if nothing else is, jobs have finished.
              echo " The number of running jobs on the queue is: "
              echo $(squeue | grep -c marcust)
              echo 
              if [[ $(squeue | grep -c marcust) -gt 1 ]]; then                                     
                  sleep 1 
              else                  
                  running=2
                  echo "running is: "
                  echo $running
                  echo
                      # Compress repeats at last point into .tgz
                      fname="data_round_"
                      round=$(echo $line4 | grep -o -P '(?<=iter).*(?=_grid)') # The round is the same whether or not the point changes.
                      pname="_point_"
                      fextension=".tar.gz"
                      Point=$(echo $line4 | grep -o -P '(?<=grid).*(?=_)')
                      Point=$(echo $Point | grep -o -P '(?<=).*(?=_)')
                      concat_name="$fname$round$pname$Point$fextension"                      

                      currDir=$(pwd)
                      cd $homeWorkFolder$dessaFolder 
                      tar czf $homeWorkFolder$dessaFolder$concat_name *.dat
                      rm *.dat # This should remove .dat files only from within this folder.
                      cd $currDir 
                      # TO Extract (within say, compObj.m), the code the same, except 'c' is replaced with 'x' for extract.
                      # currDir=$(pwd)
                      # cd $homeWorkFolder$dessaFolder 
                      # tar xzf $homeWorkFolder$dessaFolder$concat_name *.dat
                      # cd $currDir                  
              fi
          done

          sed -i '/running/d' $homeWorkFolder$jobsRunningFileOnHome             # delete 'running' variable from file
          echo running=2 >> $homeWorkFolder$jobsRunningFileOnHome             # alerts [PI_script] that all jobs have finished running
          [ -f $homeWorkFolder$jobsRunningFileOnHome ] && rm $homeWorkFolder$jobListFileOnHome        #remove [job_list.txt] from /home/ as all files have been submitted. 
          echo "job_list.txt has been deleted (in dessaJobSubmissions_Shell_Script.sh)"
          echo
          # JOBS  FINISHED RUNNING ******************************************
          #                                              ******************************************
    else
        # If I'm restarting PI_script.sh and I set running=2 manually, detect this here.
        if [ -f $homeWorkFolder$jobsRunningFileOnHome ]; then source $homeWorkFolder$jobsRunningFileOnHome; else
            echo "job_list.txt doesn't exist yet."
            echo
        fi
        echo "dessaJobSubmissions_Shell_Script.sh done for now. Sleep 4 minutes."
        echo
        sleep 240
    fi
  else
    echo "All Done!"
    echo
  fi

  else
    echo "Still waiting on resub.txt to be created."
    echo
    sleep $sleepTimeQueue    # In this case, resub.txt does not yet exist. So wait a bit.
  fi

done 

 

exit