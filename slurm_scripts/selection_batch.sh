#!/bin/bash

#SBATCH --job-name=module0_crs_sel
#SBATCH --output=output-$1-%j.txt
#SBATCH --error=output-$1-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1g
#
#SBATCH --time=2:00:00

SCRIPTS_PATH=/global/project/projectdirsdune/users/ddouglas/M0

INFILE=$1
OUTFILE=$3

COMMAND="/usr/bin/python3 $SCRIPTS_PATH/selection.py -i $INFILE -o $OUTFILE"

