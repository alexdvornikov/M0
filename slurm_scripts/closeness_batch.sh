#!/bin/bash

#SBATCH --job-name=m0_US_sel
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#
#SBATCH --constraint=haswell
#SBATCH --account=dune
#SBATCH --qos=shared

export HDF5_USE_FILE_LOCKING=FALSE

oPYTHON_EXEC=/global/common/software/nersc/cori-2022q1/sw/python/3.9-anaconda-2021.11/bin/python3
SCRIPTS_PATH=/global/project/projectdirs/dune/users/ddouglas/M0
GEOMETRY_FILE=$SCRIPTS_PATH/pixel_layouts/multi_tile_layout-2.3.16.yaml
DETECTOR_FILE=$SCRIPTS_PATH/detector_properties/module0.yaml

INFILE1=$1
INFILE2=$2
OUTFILE=$3

COMMAND="$PYTHON_EXEC $SCRIPTS_PATH/track_closeness.py $INFILE1 $INFILE2 -o $OUTFILE -d $DETECTOR_FILE -g $GEOMETRY_FILE"

$COMMAND
