SHELL:=/bin/tcsh

all:
	. /usr/share/Modules/init/tcsh; \
	module purge; \
	module load intel_comp; \
	module load gsl; \
	module load openmpi; \
	module load hdf5; \
	mpicc -o gridder_div.x gridder_div.c -lhdf5 -lm -lgmp


clean:
	rm gridder_div.x
