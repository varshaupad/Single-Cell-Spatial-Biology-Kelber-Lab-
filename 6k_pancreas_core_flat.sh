#!/bin/bash
#PBS -N scanpy_job
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -M varsha_upadhyayulla1@baylor.edu
#PBS -m abe
#PBS -o scanpy_output.log
#PBS -e scanpy_error.log

cd $PBS_O_WORKDIR

# Initialize conda in bash shell
source ~/miniconda3/etc/profile.d/conda.sh

# Activate conda environment
conda activate scanpy_env

# Run your script
SCRIPT=~/6k_pancreas_core_flat.py
python "$SCRIPT"
