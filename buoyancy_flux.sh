#!/bin/bash
#PBS -q normalbw
#PBS -P e14
#PBS -l ncpus=1
#PBS -l mem=192GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/e14+scratch/e14
#PBS -l wd
#PBS -M z5131023@unsw.edu.au

cd /g/data/e14/jb2381/CabbelingExperiments/outputs_equaldiffusion/isothermal_stepchange_nothing_660min/

# Julia
export JULIA_DEPOT_PATH="/home/561/jb2381/.julia"
export JULIA_NUM_THREADS=auto
module load julia

# Run the experiment
julia --project buoyancy_flux.jl > $PBS_JOBID.log