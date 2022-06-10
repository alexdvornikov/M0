#!/bin/bash

#SBATCH --job-name=m0_US_sel
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --constraint=haswell
#SBATCH --account=dune
#SBATCH --qos=shared

export HDF5_USE_FILE_LOCKING=FALSE

PYTHON_EXEC=/global/common/software/nersc/cori-2022q1/sw/python/3.9-anaconda-2021.11/bin/python3
SCRIPTS_PATH=/global/project/projectdirs/dune/users/olexiy/M0
GEOMETRY_FILE=$SCRIPTS_PATH/pixel_layouts/multi_tile_layout-2.3.16.yaml
DETECTOR_FILE=$SCRIPTS_PATH/detector_properties/module0.yaml

INFILE=$1
OUTFILE=$2
SELECTION_STRING=$3

COMMAND="$PYTHON_EXEC $SCRIPTS_PATH/selection.py $INFILE -o $OUTFILE -d $DETECTOR_FILE -g $GEOMETRY_FILE -c $SELECTION_STRING"

$COMMAND
