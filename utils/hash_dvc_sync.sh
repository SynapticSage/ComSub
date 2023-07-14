#!/bin/bash
# This script links the matfiles from the hash directory to the code directory
# and adds them to dvc if they are not already there

# Specify MATLAB command to be run
MATLAB_CMD="hfolder=hashdefine();disp(hfolder);"

# Run MATLAB with specified command
folder=$(matlab -nodisplay -nosplash -nodesktop -r "try, ${MATLAB_CMD}; catch, exit, end; exit")
codefolder=/Volumes/MATLAB-Drive/Shared/hash/

echo "Linking files from $folder to /Volumes/MATLAB-Drive/Shared/hash/"
echo "$folder/*.mat"
ln -sf $folder/*.mat $codefolder

for file in /Volumes/MATLAB-Drive/Shared/hash/*.mat; do
	if [ ! -f "$codefolder/${file}.dvc" ]; then
	    echo "Adding $file to dvc"
	    dvc add $codefolder/${file}
    	fi
done
