#!/bin/ksh

# This script is intended to run every Saturday at 5:00 am
#
# It updates graphs for rainfallph, rainmeter and temp_fumaroles by calling the script matlab-fonction_weekly
#
# The log directory for this script is located in: /var/log/webobs/matlab-weekly/
#
# The output of the matlab function executed by matlab-fonction_weekly is located in: /var/log/webobs/matlab-weekly/matlab_function
#
# 11/06/2022
# Arvid Ramdeane, MVO

# Setup variables, directories, filenames etc.

cd "$(dirname "$0")"

# Date and time the script is executed
DATE=$(date +%Y-%m-%d_%H-%M-%S)

# Program name of script being executed
PRGM="MATLAB-WEEKLY"
# Lock file
LOCK="/tmp/.$PRGM-lock"

# Logfile/directory setup for calling script (matlab-weekly)
# logfile variables
logfile_name="$PRGM-$DATE.log"
logfile_dir="/var/log/webobs/matlab-weekly"
# Create logfile directory if it does not exist
if [ ! -d $logfile_dir ]
then
        mkdir $logfile_dir
fi
logfile="$logfile_dir/$logfile_name"

# log directory for logging output of matlab function
matlab_log="$logfile_dir/matlab_function"
#creat directory if it does not exists
if [ ! -d $matlab_log ]
then
	mkdir $matlab_log
fi

# Directories for copying generated graphs
src_dir="/mvo/acqui"
dest_dir="/mnt/mvofls3/DVG_Data/WEBOBS_DATA"


# Begin execution of matlab functions
echo "`date`: Beginning Matlab Weekly Script." 2>&1 | tee $logfile

echo " " 2>&1 | tee -a $logfile
if [ ! -e $LOCK ]
then
        echo "`date`: Previous $LOCK does not exist" 2>&1 | tee -a $logfile
        echo " " 2>&1 | tee -a $logfile
        touch $LOCK
        #./debug_script ./matlab-fonction ovsg_pre

       	for fonction in "rainfallph(0,'all')" "rainfallph(0)" "rainmeter(0,'all')" "rainmeter(0)" "temp_fumaroles(0,'all')" "temp_fumaroles(0)"
        #for fonction in "rainfallph(0,'all')"
	do
                # Call matlab-fonction_temp to execute matlab scripts
                echo "`date`: Executing Matlab function: $fonction" 2>&1 | tee -a $logfile
		echo "`date`: Matlab output for $fonction located in: $matlab_log" 2>&1 | tee -a $logfile
                OUT="$matlab_log/$fonction.log"
		./matlab-fonction_weekly "$fonction" > $OUT

                echo " " 2>&1 | tee -a $logfile
                echo " " 2>&1 | tee -a $logfile
        done
	
        echo " " 2>&1 | tee -a $logfile

        echo "`date`: Copying graphs from $src_dir to $dest_dir" 2>&1 | tee -a $logfile

        # Check if destination directory is mounted
        MOUNTCMD="mountpoint -d $dest_dir"
        eval $MOUNTCMD
        if [ $? -eq 0 ]
        then
                #Fumarole Temps Graphs
                echo "`date`: Copying Fumarole Temps graphs." 2>&1 | tee -a $logfile
                cp "$src_dir/FumaroleTemps/Graphs/"* "$dest_dir/FumaroleTemps/Graphs/"

                # Rainfallph Graphs
                echo "`date`: Copying Rainfallph Graphs." 2>&1 | tee -a $logfile
                cp "$src_dir/RainfallpH/Graphs/"* "$dest_dir/RainfallpH/Graphs/"

                # Raingauges Graphs
                echo "`date`: Copying Raingauges Graphs" 2>&1 | tee -a $logfile
                cp "$src_dir/Raingauges/Graphs/"* "$dest_dir/Raingauges/Graphs/"

        elif [ $? -eq 1 ]
        then
                echo "`date`: Failure; incorrect invocation, permissions or system error" 2>&1 | tee -a $logfile

        else
		echo "`date`: No drives mounted at $dest_dir" 2>&1 | tee -a $logfile
        fi
        #./debug_script ./matlab-fonction ovsg_post

        rm -f $LOCK
else
        if [ $(($(date +%s)-$(find "$LOCK" -printf "%T@"|cut -d. -f1))) -lt 3600 ]
        then
                 echo "`date`: $LOCK exists. Cannot launch $PRGM, but the lock file is recent" 2>&1 | tee -a $logfile
                echo "$LOCK exists. Cannot launch $PRGM, but the lock file is recent"
        else
                echo "`date`: $LOCK exists. Cannot launch $PRGM" 2>&1 | tee -a $logfile
                echo "$LOCK exists. Cannot launch $PRGM" >&2
        fi
fi

echo "`date`: Exiting Matlab Weekly Script...." 2>&1 | tee -a $logfile


