/*
 * Parallel density gridder
 *
 * Takes 3D particle data and assigns to a uniform density grid.
 *
 * This version also calculates the divergence along each co-ordinate axis
 *
 */


#include "hdf5.h"
#include <mpi.h>
#include <gmp.h>

#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>

#include <stdlib.h>
#include "read_config_div.h"

// TODO: set as runtime value: if only this were python:-!
// PartType1: Dark Matter
#define DATASETNAME_POS "/PartType1/Coordinates"  // dataset within hdf5 file 
#define DATASETNAME_VEL "/PartType1/Velocities"  // dataset within hdf5 file 


// function initialisers
int count_files(const char *, const char *);
const char *get_filename_ext(const char *);
int offset(int, int, int, int);
//void NGP(int *, int, float, float, float);
void NGP(int *,                        // sum of number
	 double *, double *, double *,    // sum of pos
	 double *, double *, double *,    // sum of vel
	 double *, double *, double *,    // sum of pos * pos
	 double *, double *, double *,    // sum of pos * vel
	 int,                          // size of grid along each axis
	 float, float, float,          // input grid position
	 float, float, float,          // input position
	 float, float, float);         // input velocity

int main (int argc, char **argv) {

    int i,j,k;

    /*
     * Read config file
     */
    struct config_struct config;

    // read_config_file("config.txt", &config);
    read_config_file(argv[1], &config);

    char * input_directory = config.input_dir; 
    char * output_file_N = config.output_file_N;
    char * output_file_divx = config.output_file_divx;
    char * output_file_divy = config.output_file_divy;
    char * output_file_divz = config.output_file_divz;

    /*
     * Initialise weight grid
     */
    int grid_dims = config.grid_dims; //grid_size + 1;
    int slice_dim = config.slice_dim;  // width of slice in x-direction
    float sim_dims = config.sim_dims;  // simulation dimensions
    float ratio = grid_dims / sim_dims;  // ratio of grid to simulation dimensions

    hid_t 	file, dataset_pos, dataset_vel, dataspace;   // handles
    herr_t 	status;
    int		status_n;
    hsize_t 	dims[2];           			// dataset dimensions

#ifndef NO_MPI
    /*
     * Initialize MPI
     */
    int mpi_size, mpi_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    int ierr;  // store error values for MPI operations
    int root_process = 0;  // set root process to zero processor  

    if(mpi_rank == 0){
#endif
        printf("\nInput directory: %s",input_directory);
	printf("\nOutput file for density: %s", output_file_N);
	printf("\nOutput file for divx: %s", output_file_divx);
	printf("\nOutput file for divy: %s", output_file_divy);
	printf("\nOutput file for divz: %s", output_file_divz);
        printf("\nGrid dimensions: %d", grid_dims);
        printf("\nGrid slice dimension: %d", slice_dim);
        printf("\nSimulation dimensions: %lf", sim_dims);
        printf("\nConversion ratio from sim pos to grid pos: %lf\n", ratio);
#ifndef NO_MPI
    }
#endif
    /*
     * Count hdf5 files in specified directory
     */		
    const char * extension = "hdf5";
    int file_count = count_files(input_directory,extension);   
    /*
     * Find hdf5 files in specified directory
     */
    char **files = malloc(sizeof(char*) * file_count);
    int ticker = 0;
    DIR * dirp;
    struct dirent * entry;
    dirp = opendir(input_directory);
    while ((entry = readdir(dirp)) != NULL) {
	if (entry->d_type == DT_REG && !strcmp(get_filename_ext((const char *)entry->d_name), "hdf5")) {  // If the entry is a regular file, with hdf5 extension..
	    files[ticker] = calloc(sizeof(char*),sizeof(char *));  // allocate space in files array for this string
	    strcpy(files[ticker], entry->d_name);  // ...store the filename
	    strcat(files[ticker],"\0");  // add string end character to convert character array to string
	    ticker++;
	}
    }
    closedir(dirp);

#ifndef NO_MPI
    /*
     *  find number of files for given processor
     */
    int proc_files = file_count / mpi_size;
    if(mpi_rank < fmod(file_count,mpi_size)) proc_files++;
#endif
    
    /*
     * A 3D array is too big for native initialisation, and a pain using malloc.
     * So, create a 1D array and use a custom offset function (see end)
     */
    //long int grid_size = pow(grid_dims, 3); 	// Amended to calculate in slices to save memory.
    long int grid_size = grid_dims*grid_dims*slice_dim; 	
    char * fullname;
    long long int particle_count = 0;
    long long int particle_count_used = 0;
    long long int particle_count_used_total = 0;
    int *N = calloc(grid_size, sizeof *N);
    double *X = calloc(grid_size, sizeof *X);
    double *Y = calloc(grid_size, sizeof *Y);
    double *Z = calloc(grid_size, sizeof *Z);
    double *Vx = calloc(grid_size, sizeof *Vx);
    double *Vy = calloc(grid_size, sizeof *Vy);
    double *Vz = calloc(grid_size, sizeof *Vz);
    double *XX = calloc(grid_size, sizeof *XX);
    double *YY = calloc(grid_size, sizeof *YY);
    double *ZZ = calloc(grid_size, sizeof *ZZ);
    double *XVx = calloc(grid_size, sizeof *XVx);
    double *YVy = calloc(grid_size, sizeof *YVy);
    double *ZVz = calloc(grid_size, sizeof *ZVz);
    // diagnostics to check correct sim_size used
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;

    FILE *fp_N, *fp_divx, *fp_divy, *fp_divz;
#ifndef NO_MPI
    if(mpi_rank == 0){
#endif
	/*
	 * Open output files
	 */
	fp_N = fopen(output_file_N, "wb");
	if (fp_N == NULL) printf("%s could not be opened.\n",output_file_N);		
	fp_divx = fopen(output_file_divx, "wb");
	if (fp_divx == NULL) printf("%s could not be opened.\n",output_file_divx);
	fp_divy = fopen(output_file_divy, "wb");
	if (fp_divy == NULL) printf("%s could not be opened.\n",output_file_divy);
	fp_divz = fopen(output_file_divz, "wb");
	if (fp_divz == NULL) printf("%s could not be opened.\n",output_file_divz);
#ifndef NO_MPI
    }
#endif
    
    /*
     * Loop over slices
     * Ideally we would read in the data once and save, but as we are looking to minimise memory,
     * we will read in each file once for every slice of the grid.
     */
    float xgrid_min=0., xgrid_max;
    xgrid_max=(float)slice_dim;
    while (xgrid_min<0.999*grid_dims) {
	/* 
	 * Loop over iput data files
	 */
	int i_file;
#ifdef NO_MPI
	for(i_file = 0; i_file<file_count; i_file ++){
#else
	for(i_file = mpi_rank; i_file<file_count; i_file +=mpi_size){
#endif
	    //printf("\n%d %s",i_file,files[i_file]);
	    fullname = malloc(sizeof(char) * (strlen(input_directory) + strlen(files[i_file]) + 1));  // allocate space for concatenated full name and location
	    *fullname = '\0';
	    strcat(fullname, input_directory);  // concatenate directory and filename strings
	    strcat(fullname, files[i_file]);
	    //  Open the hdf5 file and dataset
	    file = H5Fopen(fullname, H5F_ACC_RDONLY, H5P_DEFAULT);
	    dataset_pos = H5Dopen(file, DATASETNAME_POS, H5P_DEFAULT);
	    dataset_vel = H5Dopen(file, DATASETNAME_VEL, H5P_DEFAULT);
	    free(fullname);
	    dataspace = H5Dget_space(dataset_pos);    // dataspace handle
	    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);  // get dataspace dimensions
	    /*
	     *  Initialise data buffer
	     */
	    int rows = dims[0];  // Number of particles?
	    int cols = dims[1];  // Number of dimensions (=3)?
	    particle_count += rows;
	    float **data_pos; 
	    float **data_vel; 
	    /* 
	     * Allocate memory for new float array[row][col] 
	     */
	    /* First allocate the memory for the top-level array (rows).
	       Make sure you use the sizeof a *pointer* to your data type. */
	    data_pos = (float**) calloc(rows, sizeof(float*));  // Done this roundabout way in case not 3 dimensional
	    data_vel = (float**) calloc(rows, sizeof(float*));
	    /* Allocate a contiguous chunk of memory for the array data values.
	       Use the sizeof the data type. */
	    data_pos[0] = (float*) calloc(cols*rows, sizeof(float));
	    data_vel[0] = (float*) calloc(cols*rows, sizeof(float));
	    /* Set the pointers in the top-level (row) array to the
	       correct memory locations in the data value chunk. */
	    for (j=1; j < rows; j++) data_pos[j] = data_pos[0] + j*cols;
	    for (j=1; j < rows; j++) data_vel[j] = data_vel[0] + j*cols;
	    /*
	     * Read dataset back.
	     */
	    status = H5Dread(dataset_pos, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data_pos[0][0]);
	    status = H5Dread(dataset_vel, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data_vel[0][0]);
	    /*
	     * Close datasets, dataspace and file
	     */
	    H5Dclose(dataset_pos);
	    H5Dclose(dataset_vel);
	    H5Sclose(dataspace);
	    H5Fclose(file);

	    /*
	     * Assign to grid.
	     * This version accumulates quantities needed to calculate the divergence
	     */
	    float xgrid, ygrid, zgrid;
	    xmin=1000.; ymin=1000.; zmin=1000.;
	    xmax=-1000.; ymax=-1000.; zmax=-1000.;
	    for(j = 0; j < rows; j++){  // loop through data rows
		xmin=fmin(xmin,data_pos[j][0]);
		xmax=fmax(xmax,data_pos[j][0]);
		ymin=fmin(ymin,data_pos[j][1]);
		ymax=fmax(ymax,data_pos[j][1]);
		zmin=fmin(zmin,data_pos[j][2]);
		zmax=fmax(zmax,data_pos[j][2]);
		xgrid = data_pos[j][0] * ratio; // x grid position
		if (xgrid<xgrid_min || xgrid>=xgrid_max) continue;
		particle_count_used++;
		xgrid -= xgrid_min;  // position relative to slice
		ygrid = data_pos[j][1] * ratio; // y grid position
		zgrid = data_pos[j][2] * ratio; // z grid position
		NGP(N,                                  // sum of number
		    X, Y, Z,                            // sum of pos
		    Vx, Vy, Vz,                         // sum of vel
		    XX, YY, ZZ,                         // sum of pos * pos
		    XVx, YVy, ZVz,                      // sum of pos * vel
		    grid_dims,                          // size of grid along each axis
		    xgrid, ygrid, zgrid,                // position in grid coordinates
		    data_pos[j][0], data_pos[j][1], data_pos[j][2],  // input position
		    data_vel[j][0], data_vel[j][1], data_vel[j][2]);  // input velocity
	    }
	    //printf("\nRank %d completed", mpi_rank);
	    free(data_pos[0]);
	    free(data_pos);
	    free(data_vel[0]);
	    free(data_vel);
	}
	/* Accumulate sums over all ranks */
#ifndef NO_MPI
	if(mpi_rank ==0){
	    MPI_Reduce(MPI_IN_PLACE, &xmin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &xmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &ymin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &ymax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &zmin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &zmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &particle_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, &particle_count_used, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, N, grid_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, X, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, Y, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, Z, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, Vx, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, Vy, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, Vz, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, XX, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, YY, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, ZZ, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, XVx, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, YVy, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(MPI_IN_PLACE, ZVz, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	} else {
	    MPI_Reduce(&xmin, &xmin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&xmax, &xmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&ymin, &ymin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&ymax, &ymax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&zmin, &zmin, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&zmax, &zmax, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&particle_count, &particle_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(&particle_count_used, &particle_count_used, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(N, N, grid_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(X, X, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(Y, Y, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(Z, Z, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(Vx, Vx, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(Vy, Vy, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(Vz, Vz, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(XX, XX, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(YY, YY, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(ZZ, ZZ, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(XVx, XVx, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(YVy, YVy, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce(ZVz, ZVz, grid_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}	    
	/*
	 * Now on the master node, we can calculate our statistics and write out
	 */
	if(mpi_rank == 0){
#endif
	    if (xgrid_min<0.001) {
		printf("x-position range of all particles is (%f, %f)\n",xmin,xmax);
		printf("y-position range of all particles is (%f, %f)\n",ymin,ymax);
		printf("z-position range of all particles is (%f, %f)\n",zmin,zmax);
	    }
	    float* divx = calloc(grid_size, sizeof *divx);
	    float* divy = calloc(grid_size, sizeof *divy);
	    float* divz = calloc(grid_size, sizeof *divz);
	    printf("Total particles: %lld\n", particle_count);
	    particle_count_used_total+=particle_count_used;
	    printf("Total particles currently processed: %lld\n", particle_count_used_total);
	    /* Calculate divergences: note that these need correcting for Hubble expansion */
	    /* ****Do I need to use doubles here: difference of two large numbers?*** */
	    for (i=0; i<grid_size; i++){
		divx[i] = (float)((N[i]*XVx[i]-X[i]*Vx[i])/(N[i]*XX[i]-X[i]*X[i]));
		divy[i] = (float)((N[i]*YVy[i]-Y[i]*Vy[i])/(N[i]*YY[i]-Y[i]*Y[i]));
		divz[i] = (float)((N[i]*ZVz[i]-Z[i]*Vz[i])/(N[i]*ZZ[i]-Z[i]*Z[i]));
	    }
	    /*
	     * Write to output files
	     */
	    fwrite(N, sizeof(int), grid_size, fp_N);
	    fwrite(divx, sizeof(float), grid_size, fp_divx);
	    fwrite(divy, sizeof(float), grid_size, fp_divy);
	    fwrite(divz, sizeof(float), grid_size, fp_divz);
	    
	    printf("Finished processing slice %f to %f.\n",xgrid_min,xgrid_max);
#ifndef NO_MPI
	}
#endif
	// Reset arrays to zero.
	particle_count=0; // Reset each slice
	particle_count_used=0;
	for (i=0; i<grid_size; i++) {
	    N[i]=0;
	    X[i]=(double)0.;
	    Y[i]=(double)0.;
	    Z[i]=(double)0.;
	    Vx[i]=(double)0.;
	    Vy[i]=(double)0.;
	    Vz[i]=(double)0.;
	    XX[i]=(double)0.;
	    YY[i]=(double)0.;
	    ZZ[i]=(double)0.;
	    XVx[i]=(double)0.;
	    YVy[i]=(double)0.;
	    ZVz[i]=(double)0.;
	}
	/*
	 * Update slice counters
	 */
	xgrid_min+=(float)slice_dim;
	xgrid_max=xgrid_min+(float)slice_dim;
	// End of loop over slices
    }
    // Close output files
#ifndef NO_MPI
    if(mpi_rank == 0){
#endif
	fclose(fp_N);
	fclose(fp_divx);
	fclose(fp_divy);
	fclose(fp_divz);
#ifndef NO_MPI
    }
#endif
    // Let's be good and free up all the memory, even though the end of the program.
    free(N);
    free(X);
    free(Y);
    free(Z);
    free(Vx);
    free(Vy);
    free(Vz);
    free(XX);
    free(YY);
    free(ZZ);
    free(XVx);
    free(YVy);
    free(ZVz);
    free(files[0]);
    free(files);
#ifndef NO_MPI
    ierr = MPI_Finalize();
#endif
    return 0;
}


/*
 * Nearest Grid Point assignment
 */
void NGP(int   *N,                              // sum of number
	 double *X, double *Y, double *Z,          // sum of pos
	 double *Vx, double *Vy, double *Vz,       // sum of vel
	 double *XX, double *YY, double *ZZ,       // sum of pos * pos
	 double *XVx, double *YVy, double *ZVz,    // sum of pos * vel
	 int dims,                              // size of grid along each axis
	 float xgrid, float ygrid, float zgrid, // input grid position
	 float x, float y, float z,             // input position
	 float vx, float vy, float vz){         // input velocity
    int loc;
    float xrel, yrel, zrel;
    //int loc_trace=1118; // Used to print diagnostic info to debug
    loc = offset((int) xgrid, (int) ygrid, (int) zgrid, dims);   // location in 1-d array
    //if (loc==loc_trace) printf("loc %d: pos (%f, %f, %f), vel (%f, %f, %f)\n",loc,x,y,z,vx,vy,vz);
    N[loc] += 1;
    // Measure relative positions from 0 to 1 to avoid large numbers in the accumulating sums below
    xrel = x-(int)xgrid;
    yrel = y-(int)ygrid;
    zrel = z-(int)zgrid;
    X[loc] += xrel;
    Y[loc] += yrel;
    Z[loc] += zrel;
    Vx[loc] += vx;
    Vy[loc] += vy;
    Vz[loc] += vz;
    XX[loc] += xrel*xrel;
    YY[loc] += yrel*yrel;
    ZZ[loc] += zrel*zrel;
    XVx[loc] += xrel*vx;
    YVy[loc] += yrel*vy;
    ZVz[loc] += zrel*vz;
}

/*
 * Cloud In Cell assignment
 */
void CIC(int * array, int dims, int x, int y, int z){
    
    array[offset((int)x, (int)y, (int)z, dims)] += (1 - fmod(x, 1.)) * (1 - fmod(y, 1.)) * (1 - fmod(z, 1.));
    array[offset((int)x + 1, (int)y, (int)z, dims)] += fmod(x, 1.) * (1 - fmod(y, 1.)) * (1 - fmod(z, 1.));
    
    array[offset((int)x, (int)y + 1, (int)z, dims)] += (1 - fmod(x, 1.)) * fmod(y, 1.) * (1 - fmod(z, 1.));
    array[offset((int)x + 1, (int)y + 1, (int)z, dims)] += fmod(x, 1.) * fmod(y, 1.) * (1 - fmod(z, 1.));
    
    array[offset((int)x, (int)y, (int)z + 1, dims)] += (1 - fmod(x, 1.)) * (1 - fmod(y, 1.)) * fmod(z, 1.);
    array[offset((int)x + 1, (int)y, (int)z + 1, dims)] += fmod(x, 1.) * (1 - fmod(y, 1.)) * fmod(z, 1.);
     
    array[offset((int)x, (int)y + 1, (int)z + 1, dims)] += (1 - fmod(x, 1.)) * fmod(y, 1.) * fmod(z, 1.);
    array[offset((int)x + 1, (int)y + 1, (int)z + 1, dims)] += fmod(x, 1.) * fmod(y, 1.) * fmod(z, 1.);  
 
}

/*
 * Triangular Shaped Cloud assignment
 */
void TSC(int * array, int dims, int x, int y, int z){
}


int count_files(const char * directory, const char * extension){
	/*
	 * given a directory, count the number of files in it
	 */

	int file_count = 0;
	DIR * dirp = opendir(directory);
	struct dirent * entry;

	while ((entry = readdir(dirp)) != NULL) {
		if (entry->d_type == DT_REG && !strcmp(get_filename_ext((const char *)entry->d_name), extension)) file_count++;  // If the entry is a regular file..
	}
	closedir(dirp);

	return file_count;
}


const char *get_filename_ext(const char * filename) {
	/*
	 * Given a filename, return the extension
	 */

    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}


/*
 * Given a 3D array (grid_dims^3) flattened to 1D, return offset for given 3D coordinates in flat array
 */
/* 
 * Incorrect F-style ordering
 int offset(int x, int y, int z, int grid_dims) { return ( z * grid_dims * grid_dims ) + ( y * grid_dims ) + x ; }
*/
int offset(int x, int y, int z, int grid_dims) {
    return (x * grid_dims + y) * grid_dims + z ; }
