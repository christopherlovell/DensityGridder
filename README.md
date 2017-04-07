
# Density Gridder

Takes 3D particle data and maps on to a density grid, using MPI. Can optionally apply convolution using serial R scripts.

Given a box size `lbox`, and grid size `grid`, the gridder finds the density using a Nearest Grid Point assignment scheme.

|---------`lbox`----------|
|-|-|-|-|-|-|-|-|-|-|-|-| <- `grid`=12

Center of each grid point in each dimension = `lbox/grid/2`

## Configuration

Configuration parameters are contained in config.txt, which now does not need to be set at compile time. The configurable parameters are as follows:

INPUT_DIRECTORY: location of input particle files
OUTPUT_DIRECTORY: location of output file, with name and extension
GRID_DIMS: grid dimensions over which density calculated (`grid`)
SIM_DIMS: size of simulation box along one side (`lbox`)


## Compile

To compile gridder.c execute the following in the top level directory.

```
make clean
make
```

## To run 

To run gridder on a Platform LSF batch system, submit grid_run.sh. You may need to change the modules depending on your system availability / compiler requirements. To run on the login node you need to load the following modules.

```
module unload intel_comp
module load gnu_comp
module load openmpi
module load hdf5

mpirun -np 4 gridder.x
```
