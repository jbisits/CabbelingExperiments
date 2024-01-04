#!/bin/bash
#PBS -q gpuvolta
#PBS -P e14
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l jobfs=1GB
#PBS -l walltime=32:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/CabbelingExperiments/

# Julia
export JULIA_DEPOT_PATH="/home/561/jb2381/.julia"
export JULIA_NUM_THREADS=auto
module load julia

# Run the experiment
# julia --project cabbeling_gpu_salinitynoise.jl > $PBS_JOBID.log
julia --project cabbeling_gpu_salinitynoise.jl > $PBS_JOBID.log