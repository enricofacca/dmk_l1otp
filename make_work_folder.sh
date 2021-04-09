#!/bin/bash

# store repository location
pathrepo=$(readlink -e $(pwd)/../)

###############################################
# Make working directory with controls and runs
###############################################
fname=$1
if [ -z "$fname" ]
then
    echo 'Give path of working folder'
    exit 1
fi
# create working folder
mkdir ${fname} 
# move 
mkdir ${fname}/runs
echo ${pathrepo}>${fname}/location_repository.txt
cp -r ${pathrepo}/tools/script_controls/dmk.inputs     ${fname}
cp -r ${pathrepo}/tools/script_controls/dmk.ctrl       ${fname}
cp -r ${pathrepo}/tools/script_controls/dmk.inputs     ${fname}
cp -r ${pathrepo}/tools/script_controls/dmk_folder.py  ${fname}
