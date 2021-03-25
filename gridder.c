/*
 * Parallel density gridder
 *
 * Takes 3D particle data and assigns to a uniform density grid.
 */


#include "hdf5.h"
#include <mpi.h>
#include <gmp.h>

#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>

#include <stdlib.h>
#include "read_config.h"

// TODO: set as runtime value
// PartType1: Dark Matter
#define DATASETNAME "/PartType1/Coordinates"  // dataset within hdf5 file 


// function initialisers
int count_files(const char *, const char *);
const char *get_filename_ext(const char *);
int offset(int, int, int, int);
void NGP(int *, int, float, float, float);

int main (int argc, char **argv) {

    int i,j,k;

    /*
     * Read config file
     */
    struct config_struct config;

    // read_config_file("config.txt", &config);
    read_config_file(argv[1], &config);

	char * input_directory = config.input_dir; 
	char * output_file = config.output_dir;

	/*
	 * Initialise weight grid
	 */
	int grid_dims = config.grid_dims; //grid_size + 1; 
	double sim_dims = config.sim_dims;  // simulation dimensions

	hid_t 		file, dataset, dataspace;   // handles
	herr_t 		status;
	int		status_n;
	hsize_t 	dims[2];           			// dataset dimensions
	
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
        printf("\nInput directory: %s\nOutput directory+filename: %s", input_directory, output_file);
        printf("\nGrid dimensions: %d", grid_dims);
        printf("\nSimulation dimensions: %lf", sim_dims);
    }
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


	/*
	 *  find number of files for given processor
	 */
	int proc_files = file_count / mpi_size;
	if(mpi_rank < fmod(file_count,mpi_size)) proc_files++;


	/*
	 * A 3D array is too big for native initialisation, and a pain using malloc.
	 * So, create a 1D array and use a custom offset function (see end)
	 */
	long int w_grid_size = pow(grid_dims, 3); 
	
	int *w = calloc(w_grid_size, sizeof *w);
	int *w_slave = calloc(w_grid_size, sizeof *w_slave);
	
	char * fullname;

    long long int particle_count = 0;
    long long int particle_count_slave = 0;

	for(i = mpi_rank; i<file_count; i+=mpi_size){

		printf("%d %s\n",i,files[i]);

		fullname = malloc(sizeof(char) * (strlen(input_directory) + strlen(files[i]) + 1));  // allocate space for concatenated full name and location
		*fullname = '\0';

		strcat(fullname, input_directory);  // concatenate directory and filename strings
		strcat(fullname, files[i]);

		//  Open the hdf5 file and dataset
		file = H5Fopen(fullname, H5F_ACC_RDONLY, H5P_DEFAULT);
		dataset = H5Dopen(file, DATASETNAME, H5P_DEFAULT);

		free(fullname);

		dataspace = H5Dget_space(dataset);    // dataspace handle
		status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);  // get dataspace dimensions

		/*
		 *  Initialise data buffer
		 */
		int rows = dims[0];
		int cols = dims[1];

        particle_count_slave += rows;

        printf("%d: %d particles, %lld total\n", i, rows, particle_count_slave);
		
		float **data_out; 
		
		/* 
		 * Allocate memory for new float array[row][col] 
		 */

		/* First allocate the memory for the top-level array (rows).
		Make sure you use the sizeof a *pointer* to your data type. */
		data_out = (float**) calloc(rows, sizeof(float*));

		/* Allocate a contiguous chunk of memory for the array data values.
		Use the sizeof the data type. */
		data_out[0] = (float*) calloc(cols*rows, sizeof(float));

		/* Set the pointers in the top-level (row) array to the
		correct memory locations in the data value chunk. */
		for (j=1; j < rows; j++) data_out[j] = data_out[0] + j*cols;

		/*
		 * Read dataset back.
		 */
		status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data_out[0][0]);


		H5Dclose(dataset);
		H5Sclose(dataspace);
		H5Fclose(file);

		/*
		 * Assign to grid
		 */
		float xpos, ypos, zpos;

		double ratio = grid_dims / sim_dims;  // ratio of grid to simulation dimensions
		
		for(j = 0; j < rows; j++){  // loop through data rows 
			
			xpos = data_out[j][0] * ratio; // x grid position
			ypos = data_out[j][1] * ratio; // y grid position
			zpos = data_out[j][2] * ratio; // z grid position

			// w_slave[offset((int)xpos, (int)ypos, (int)zpos, grid_dims)] += 1;
            
            NGP(w_slave, grid_dims, xpos, ypos, zpos);  // NGP assignment (CIC & TSC to be developed)

		}

        printf("%d complete\n", i);
			
		free(data_out[0]);
		free(data_out);
	}
	
	MPI_Reduce(w_slave, w, w_grid_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&particle_count_slave, &particle_count, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	
    if(mpi_rank == 0){

        printf("Total particles: %lld\n", particle_count);

        /*
         * Check particle number
         */
        /*
        mpz_t sum;
        mpz_t temp;

        mpz_init(sum);
        mpz_init(temp);

        gmp_printf("%Zd %Zd",sum,temp);

        for(i=0; i < w_grid_size; i++){
            if(i % 1000000 == 0){
                gmp_printf("%Zd \n",sum);
            }                        

            mpz_set_si(temp,w[i]);
            mpz_add(sum,sum,temp);
        }

        gmp_printf("%Zd %d %d %d\n",sum,w[w_grid_size-1],w[w_grid_size],w[0]);

        mpz_clear(sum);
        mpz_clear(temp);

        */

        /*
         * For periodic box, combine grid edges 
         */
/*
        printf("Combining grid edges.\n");

        for(i = 0; i < grid_dims; i++){
            for(j = 0; j < grid_dims; j++){
                w[offset(0, i, j, grid_dims)] += 1; 
                //w[offset(grid_size, i, j, grid_dims)];
                w[offset(i, 0, j, grid_dims)] += 1;
                //w[offset(i, grid_size, j, grid_dims)];
                w[offset(i, j, 0, grid_dims)] += 1;
                //w[offset(i, j, grid_size, grid_dims)];
            }
        }
*/

        /*
         * Write to file
         */

		FILE *fp;		
		fp = fopen(output_file, "wb");
		
		if (fp == NULL) printf("File could not be opened.\n");		
	
		// Save final weight array
		for(i = 0; i < grid_dims; i++){
			for(j = 0; j < grid_dims; j++){
				for(k = 0; k < grid_dims; k++){
					fprintf(fp, "%d\n",  w[offset(k, j, i, grid_dims)]);
				}

			    //fprintf(fp, "\n");
			}
		    //fprintf(fp, "\n");
		}
	    
        fclose(fp);

	}

	free(w);
	free(w_slave);

	free(files[0]);
	free(files);

	ierr = MPI_Finalize();

	return 0;
}


/*
 * Nearest Grid Point assignment
 */
void NGP(int * array, int dims, float x, float y, float z){

    array[offset((int) x, (int) y, (int) z, dims)] += 1;

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
int offset(int x, int y, int z, int grid_dims) { return ( z * grid_dims * grid_dims ) + ( y * grid_dims ) + x ; }



