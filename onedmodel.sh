#!/bin/bash
#PBS -q normalbw
#PBS -P e14
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/CabbelingExperiments/1DModel

# Julia
export JULIA_DEPOT_PATH="/home/561/jb2381/.julia"
export JULIA_NUM_THREADS=auto
module load julia

# Run the experiment
julia --project run_onedmodel.jl > $PBS_JOBID.log