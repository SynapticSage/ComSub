#!/bin/bash
echo "Backing up shared MATLAB drive to citadel"
rsync -avu /Volumes/MATLAB-Drive/Shared citadel:/volume1/sharespace-commsub/
echo "Backing up shared dvc drive to citadel"
rsync -avu /Volumes/RY20-RY9-RY7-dlc-drive/commsubspace_dvc citadel:/volume1/sharespace-commsub/
