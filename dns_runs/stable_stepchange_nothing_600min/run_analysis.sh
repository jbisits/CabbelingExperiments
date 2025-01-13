#!/bin/bash
#PBS -q normal
#PBS -P e14
#PBS -l ncpus=1
#PBS -l mem=128GB
#PBS -l walltime=08:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/CabbelingExperiments/dns_runs/stable_stepchange_nothing_600min

# Julia
export JULIA_DEPOT_PATH="/g/data/e14/jb2381/.julia"
export JULIA_NUM_THREADS=auto

# Run the experiment
julia --project analysis_and_saving.jl > $PBS_JOBID.log