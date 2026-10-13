#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CONFIG_VARIABLE_LEN 1024
#define CONFIG_LINE_BUFFER_SIZE 100
#define MAX_LLIST_NAME_LEN 1024
#define MAX_OUT_NAME_LEN 1024

//struct config_struct
//typedef struct 
struct config_struct { 
    char input_dir[MAX_LLIST_NAME_LEN];
    char output_file_N[MAX_LLIST_NAME_LEN];
    char output_file_divx[MAX_LLIST_NAME_LEN];
    char output_file_divy[MAX_LLIST_NAME_LEN];
    char output_file_divz[MAX_LLIST_NAME_LEN];
    int grid_dims;
    int slice_dim;
    double sim_dims;
};

void read_int_from_config_line(char* config_line, int* val) { 
    char prm_name[MAX_CONFIG_VARIABLE_LEN];
    sscanf(config_line, "%s %i\n", prm_name, val);
}

void read_double_from_config_line(char* config_line, double* val) {    
    char prm_name[MAX_CONFIG_VARIABLE_LEN];
    sscanf(config_line, "%s %lf\n", prm_name, val);
}

void read_str_from_config_line(char* config_line, char* val) {    
    char prm_name[MAX_CONFIG_VARIABLE_LEN];
    sscanf(config_line, "%s %s\n", prm_name, val);
}


void read_config_file(char* config_filename, struct config_struct* conf) {

    FILE *fp;
    char buf[CONFIG_LINE_BUFFER_SIZE];

    if ((fp=fopen(config_filename, "r")) == NULL) {
        fprintf(stderr, "Failed to open config file %s", config_filename);
        exit(EXIT_FAILURE);
    }

    while(! feof(fp)) {
        fgets(buf, CONFIG_LINE_BUFFER_SIZE, fp);
        if (buf[0] == '#' || strlen(buf) < 4) {
            continue;
        }
        if (strstr(buf, "INPUT_DIR ")) {
            read_str_from_config_line(buf, conf->input_dir);
        }
        if (strstr(buf, "OUTPUT_FILE_N ")) {
            read_str_from_config_line(buf, conf->output_file_N);
        }
        if (strstr(buf, "OUTPUT_FILE_DIVX ")) {
            read_str_from_config_line(buf, conf->output_file_divx);
        }
        if (strstr(buf, "OUTPUT_FILE_DIVY ")) {
            read_str_from_config_line(buf, conf->output_file_divy);
        }
        if (strstr(buf, "OUTPUT_FILE_DIVZ ")) {
            read_str_from_config_line(buf, conf->output_file_divz);
        }
        if (strstr(buf, "GRID_DIMS ")) {
            read_int_from_config_line(buf, &conf->grid_dims);
        }
        if (strstr(buf, "SLICE_DIM ")) {
            read_int_from_config_line(buf, &conf->slice_dim);
        }
        if (strstr(buf, "SIM_DIMS ")) {
            read_double_from_config_line(buf, &conf->sim_dims);
        }
    }

    fclose(fp);
}

