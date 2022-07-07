#!/bin/bash

#SBATCH --job-name=m0_US_sel
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#SBATCH --time=10:00:00
#SBATCH --constraint=haswell
#SBATCH --account=dune
#SBATCH --qos=shared

export HDF5_USE_FILE_LOCKING=FALSE

conda activate /global/common/software/dune/module0_flow_nompi

SCRIPTS_PATH=/global/project/projectdirs/dune/users/olexiy/M0
GEOMETRY_FILE=$SCRIPTS_PATH/pixel_layouts/module1_layout-2.3.16.yaml
DETECTOR_FILE=$SCRIPTS_PATH/detector_properties/module0.yaml

INFILE=$1
OUTFILE1=$2
OUTFILE2=$3
OUTFILE3=$4

# COMMAND="python $SCRIPTS_PATH/selection_m1.py $INFILE -o $OUTFILE -d $DETECTOR_FILE -g $GEOMETRY_FILE"
COMMAND="python $SCRIPTS_PATH/selection_curvature.py $INFILE -o1 $OUTFILE1 -o2 $OUTFILE2 -o3 $OUTFILE3 -d $DETECTOR_FILE -g $GEOMETRY_FILE"

$COMMAND