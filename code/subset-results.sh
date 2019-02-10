#!/bin/bash

## ################################################################
## This script takes a file name as the first argument {STRING}
## and the file path as a second argument {DIR_SCRATCH} and outputs
## only rows that have a p-value < 1e-8. P-values have to be in the
## 5th column.

## The script can be called through R using the following command
## system(glue::glue("./subset-results.sh {STRING} {DIR_SCRATCH}"))
## ################################################################

STRING=$1 ## takes first argument
DIR_SCRATCH=$2 ## takes first argument

## select all model output files
FILES=$(ls $DIR_SCRATCH/lmm_gaston_*$STRING* | grep -v 'subset')

## this loops through all files, creates a new file name and
## only writes out rows with the 5th column having a value < 1e-8.
for FILENAME in $FILES
do
  NEWFILENAME="$FILENAME-subset.txt" 
  
  if [ ! -f $NEWFILENAME ]; then
    echo $NEWFILENAME
    awk '{ if ($5 < 1e-8) { print } }' $FILENAME > $NEWFILENAME &
  fi
done
