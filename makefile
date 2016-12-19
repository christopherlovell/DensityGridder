SHELL:=/bin/tcsh

all:
	. /usr/share/Modules/init/tcsh; \
    module purge; \
	module load intel_comp/c4/2013.1.046 gsl/1.16 platform_mpi/9.1.2 hdf5/1.8.12; \
	mpicc -I/usr/include/hdf5/serial -o gridder.x -L/usr/lib/x86_64-linux-gnu/hdf5/serial gridder.c -lhdf5 -lm -lgmp

clean:
	rm gridder.x
