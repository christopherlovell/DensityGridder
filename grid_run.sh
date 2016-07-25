#!/bin/tcsh
#BSUB -L /bin/tcsh
#BSUB -q bench1
#BSUB -P durham
#BSUB -J density_gridder
#BSUB -W 00:30
#BSUB -n 12
#BSUB -o logs/grid_out.%J
#BSUB -e logs/grid_err.%J

module purge
module load intel_comp/c4/2013.1.046 gsl/1.16 platform_mpi/9.1.2 hdf5/1.8.12

mpirun gridder.x

