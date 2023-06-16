#!/bin/bash
# analysis.sh - an analysis program
# $1 (input) and $2 (output) are the first and second arguments to this script

# strip off the directory paths to get just the filename
# BASE=`basename $1`
varname=$1
sitename=$2
seed=$3

echo 
echo "run('$1');"
date

# Run script
#R --no-save < $1
Rscript ./scripts/02_run_model.R $varname $sitename $seed
echo "run('$1'); done"
date
