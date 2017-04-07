SHELL:=/bin/tcsh

## Alternative compilation options:
# module load intel_comp/c4/2013.1.046 gsl/1.16 platform_mpi/9.1.2 hdf5/1.8.12; 
# mpicc -I/usr/include/hdf5/serial -o gridder.x -L/usr/lib/x86_64-linux-gnu/hdf5/serial gridder.c -lhdf5 -lm -lgmp

all:
	. /usr/share/Modules/init/tcsh; \
    module purge; \
	module unload intel_comp/c4/2015; \
	module unload gnu_comp; \
	module unload platform_mpi; \
	module load intel_comp/c5/2013.0.028; \
	module load hdf5/1.8.9; \
	module load openmpi; \
	module load hdf5; \
	mpicc -o gridder.x gridder.c -lhdf5 -lm -lgmp


clean:
	rm gridder.x
