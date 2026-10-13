#!/bin/bash
#SBATCH -p cosma7-shm
#SBATCH -A dp004
#SBATCH --job-name=density_gridder_div_eagle
#SBATCH --output=logs/grid_out_ref.%J
#SBATCH --error=logs/grid_err_ref.%J
#SBATCH -t 240:00
#SBATCH --ntasks 28
#SBATCH --exclusive

module purge
module load intel_comp
module load gsl
module load openmpi
module load hdf5 

mpicc -g -check-pointers=rw -o gridder_div.x gridder_div.c -lhdf5 -lm -lgmp

mpirun gridder_div.x configs/eagle_config.txt

