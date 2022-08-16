#!/bin/bash

#SBATCH --job-name=cath_cross_%j
#SBATCH --output=output-%j.txt
#SBATCH --error=output-%j.txt
#SBATCH --time=4:00:00
#SBATCH --constraint=haswell
#SBATCH --account=dune
#SBATCH --qos=shared

export HDF5_USE_FILE_LOCKING=FALSE

module load python
conda activate /global/common/software/dune/module0_flow_nompi

SCRIPTS_PATH=/global/project/projectdirs/dune/users/ddouglas/M0
GEOMETRY_FILE=$SCRIPTS_PATH/pixel_layouts/module1_layout-2.3.16.yaml
DETECTOR_FILE=$SCRIPTS_PATH/detector_properties/module0.yaml

INFILE=$1
OUTFILE=$2

COMMAND="python $SCRIPTS_PATH/selection_cathode_crossers_volumetric.py $INFILE -o $OUTFILE -d $DETECTOR_FILE -g $GEOMETRY_FILE"

$COMMAND
