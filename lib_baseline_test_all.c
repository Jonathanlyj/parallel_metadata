/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>
#include "baseline_ncx_lib.h"
#include <malloc.h>



static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

#define SOURCE_NAME "/global/homes/y/yll6162/parallel_metadata/data/nue_slice_panoptic_hdf_merged.nc"
// #define SOURCE_NAME "/pscratch/sd/y/yll6162/FS_2M_8/nue_slice_panoptic_hdf_merged_10_copy.nc"
// #define SOURCE_NAME "/files2/scratch/yll6162/parallel_metadata/script/dummy_test.nc"
#define OUTPUT_NAME "/pscratch/sd/y/yll6162/FS_2M_8/lib_baseline_test_all.nc"
// #define OUTPUT_NAME "/pscratch/sd/y/yll6162/FS_2M_32/lib_baseline_test_all.nc"
//  #define OUTPUT_NAME "lib_baseline_test_all.nc"

double def_start_time;
double total_def_time = 0;

/*----< pnetcdf_check_mem_usage() >------------------------------------------*/
/* check PnetCDF library internal memory usage */
static int
pnetcdf_check_mem_usage(MPI_Comm comm)
{
    int err, nerrs=0, rank;
    MPI_Offset malloc_size, sum_size;

    MPI_Comm_rank(comm, &rank);

    /* print info about PnetCDF internal malloc usage */
    err = ncmpi_inq_malloc_max_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0)
            printf("maximum heap memory allocated by PnetCDF internally is %lld bytes (%.2f MB)\n",
                   sum_size, (float)sum_size /1048576);
    }
    else if (err != NC_ENOTENABLED) {
        printf("Error at %s:%d: %s\n", __FILE__,__LINE__,ncmpi_strerror(err));
        nerrs++;
    }
    return nerrs;
}


void print_memory_info() {
    struct mallinfo info = mallinfo();
    printf("\nMemory Allocation Info:\n");
    printf("Total space in arena: %zu bytes (%.2f MB) \n", info.arena, (double)info.arena /1048576);
    printf("Total allocated space: %zu bytes (%.2f MB)\n", info.uordblks, (double)info.uordblks /1048576);
    printf("Total free space: %zu bytes (%.2f MB)\n", info.fordblks, (double)info.fordblks /1048576); 
    printf("Largest free block: %zu bytes (%.2f MB)\n", info.keepcost, (double)info.keepcost /1048576);
    printf("------------------------------------\n");
}
/* ---------------------------------- Read Metadata ----------------------------------------*/


void read_metadata(int rank, int size, struct hdr *file_info) {

    int ncid, num_vars, num_dims, tot_num_dims, elem_sz, v_attrV_xsz, status;

    MPI_Offset start, count;
    file_info->dims.ndefined = 0;
    file_info->attrs.ndefined = 0;
    file_info->vars.ndefined = 0;
    file_info->xsz = 0;
    tot_num_dims = 0;
    // Open the NetCDF file
    if (ncmpi_open(MPI_COMM_WORLD, SOURCE_NAME, NC_NOWRITE, MPI_INFO_NULL, &ncid) != NC_NOERR) {
        fprintf(stderr, "Error opening NetCDF file.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Get the number of variables
    if (ncmpi_inq_nvars(ncid, &num_vars) != NC_NOERR) {
        fprintf(stderr, "Error getting the number of variables.\n");
        ncmpi_close(ncid);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    file_info->xsz += 2 * sizeof(uint32_t); // NC_Variable and ndefined
    file_info->xsz += 2 * sizeof(uint32_t); // NC_Dimension and nelems
    // Calculate equal distribution of variables among processes
    int vars_per_process = num_vars / size;
    int remainder = num_vars % size;
    // int remainder = num_vars % 8;    
    // int vars_per_process = num_vars / 8;
    // Determine start and count based on rank
    start = rank * vars_per_process + (rank < remainder ? rank : remainder);
    count = vars_per_process + (rank < remainder ? 1 : 0);
    // if (rank == 0){
    //     printf("\nNumber of variables per process: %d to %d\n", vars_per_process, vars_per_process + 1);
    // }

    file_info->vars.ndefined = count;
    file_info->vars.value = (hdr_var **)malloc(file_info->vars.ndefined * sizeof(hdr_var *));
    // Each process reads its subset of variables

    for (int i = start; i < start + count; ++i) {
        // Get variable information
        hdr_var *variable_info = (hdr_var *)malloc(sizeof(hdr_var));
        variable_info->ndims = 0;  // Initialize the number of dimensions
        variable_info->attrs.ndefined = 0;  // Initialize the number of attributes

        // Get variable information
        char var_name[NC_MAX_NAME + 1];
        ncmpi_inq_varname(ncid, i, var_name);
        variable_info->name_len = strlen(var_name);
        variable_info->name = (char *)malloc((variable_info->name_len + 1) * sizeof(char));
        strcpy(variable_info->name, var_name);

        // Get number of dimensions
        ncmpi_inq_varndims(ncid, i, &(variable_info->ndims));

        num_dims = variable_info->ndims;

        // Allocate memory for dimension IDs
        variable_info->dimids = (int *)malloc(variable_info->ndims * sizeof(int));

        // Get dimension IDs
        ncmpi_inq_vardimid(ncid, i, variable_info->dimids);

        // Read attributes
        nc_type xtype;
        int num_attrs;
        ncmpi_inq_var(ncid, i, NULL, &xtype, NULL, NULL, &num_attrs);
        variable_info->attrs.ndefined = num_attrs;
        variable_info->xtype = xtype;

        file_info->xsz += sizeof(uint32_t) + sizeof(char) * variable_info -> name_len; //var name
        file_info->xsz += sizeof(uint32_t); //xtype
        file_info->xsz += sizeof(uint32_t); //nelems of dim list
        file_info->xsz += sizeof(uint32_t) * variable_info ->ndims; // dimid list


    
        //Read and store dimension information 
        file_info->dims.ndefined += num_dims;
        if (tot_num_dims == 0) {
            file_info->dims.value = (hdr_dim **)malloc(file_info->dims.ndefined * sizeof(hdr_dim *));
        }else{
            file_info->dims.value = (hdr_dim **)realloc(file_info->dims.value, file_info->dims.ndefined * sizeof(hdr_dim *));
        }
       
        for (int k = 0; k < num_dims; ++k) {
            hdr_dim *dimension_info = (hdr_dim *)malloc(sizeof(hdr_dim));
            int dimid = variable_info->dimids[k];
 
                
            // Get dimension name
            char dim_name[NC_MAX_NAME + 1];
            ncmpi_inq_dimname(ncid, dimid, dim_name);
            dimension_info->name_len = strlen(dim_name);
            dimension_info->name = (char *)malloc((dimension_info->name_len + 1) * sizeof(char));
            strcpy(dimension_info->name, dim_name);
            // if (i - start <3){
            //     printf("rank:%d varid: %d, dimid: %d, dimname: %s\n", rank, i, dimid, dim_name);
            // }   
            // Get dimension size
            ncmpi_inq_dimlen(ncid, dimid, &(dimension_info->size));

            
            file_info->dims.value[k + tot_num_dims] = dimension_info;
            variable_info->dimids[k] = k + tot_num_dims; //overwriting previous global dim id to local dim id
            file_info->xsz += sizeof(uint32_t) + sizeof(char) * dimension_info -> name_len; // dim name
            file_info->xsz += sizeof(uint32_t); //size

        }
        tot_num_dims += num_dims;

        // Allocate memory for attributes
        file_info->xsz += 2 * sizeof(uint32_t); // NC_Attribute and ndefine
        variable_info->attrs.value = (hdr_attr **)malloc(num_attrs * sizeof(hdr_attr *));
        for (int j = 0; j < num_attrs; ++j) {
            variable_info->attrs.value[j] = (hdr_attr *)malloc(sizeof(hdr_attr));
            variable_info->attrs.value[j]->name_len = 0;  // Initialize attribute name length

            // Get attribute name 
            char att_name[NC_MAX_NAME + 1];
            ncmpi_inq_attname(ncid, i, j, att_name);
            variable_info->attrs.value[j]->name_len = strlen(att_name);
            variable_info->attrs.value[j]->name = (char *)malloc(( strlen(att_name) + 1) * sizeof(char));
            strcpy(variable_info->attrs.value[j]->name, att_name);

            // Get attribute type and size
            nc_type attr_type;
            MPI_Offset attr_size;
            ncmpi_inq_att(ncid, i, variable_info->attrs.value[j]->name, &attr_type, &attr_size);

            // Allocate memory for attribute value and read it
            variable_info->attrs.value[j]->xtype = attr_type;
            variable_info->attrs.value[j]->nelems = attr_size;
            meta_xlen_nc_type(attr_type, &elem_sz);
            variable_info->attrs.value[j]->xvalue = malloc(attr_size * elem_sz);
            ncmpi_get_att(ncid, i, variable_info->attrs.value[j]->name, variable_info->attrs.value[j]->xvalue);

            file_info->xsz += sizeof(uint32_t) + sizeof(char) * variable_info->attrs.value[j]->name_len; //attr name
            file_info->xsz += sizeof(uint32_t); // nc_type
            file_info->xsz += sizeof(uint32_t); // nelems
            status = meta_xlen_nc_type(variable_info->attrs.value[j]->xtype, &v_attrV_xsz);
            file_info->xsz += variable_info->attrs.value[j]->nelems * v_attrV_xsz; // attr_value
        }

        // Add the variable information to the file_info structure

        file_info->vars.value[i - start] = variable_info;
    }

    // Close the source NetCDF file
    ncmpi_close(ncid);
}

/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr(struct hdr *hdr_data, int ncid, int rank){
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int *dimid = (int *)malloc(ndims * sizeof(int));
    int i,j,k,nerrs=0;
    int err;

    for (i=0; i<ndims; i++){
        def_start_time = MPI_Wtime();
        err = ncmpi_def_dim(ncid, hdr_data->dims.value[i]->name,  hdr_data->dims.value[i]->size, &dimid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        // }
    }
    //define variables
    int nvars = hdr_data->vars.ndefined;
    int *varid = (int *)malloc(nvars * sizeof(int));
    int v_ndims, v_namelen, xtype, n_att;
    int *v_dimids;
    int att_namelen, att_xtype, att_nelems;

    for (i=0; i<nvars; i++){

        v_namelen =  hdr_data->vars.value[i]->name_len;
        xtype = hdr_data->vars.value[i]->xtype;

        v_ndims = hdr_data->vars.value[i]->ndims;
        v_dimids = (int *)malloc(v_ndims * sizeof(int));
        for(j=0; j<v_ndims; j++) v_dimids[j] = dimid[hdr_data->vars.value[i]->dimids[j]];
        def_start_time = MPI_Wtime();
        err = ncmpi_def_var(ncid, hdr_data->vars.value[i]->name, xtype, v_ndims,  v_dimids, &varid[i]); ERR
        total_def_time += MPI_Wtime() - def_start_time;
        n_att = hdr_data->vars.value[i]->attrs.ndefined;
        // printf("\nn_att: %d\n", n_att);
        for(k=0; k<n_att; k++){
            att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
            att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
            att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
            err = ncmpi_put_att(ncid, varid[i],  hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
            att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
        }

    }
    return nerrs;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size, status, err, nerrs=0;
    double end_time, start_time, start_time1, end_time1, start_time2, start_time3, max_time, min_time;
    double meta_create_time, enddef_time, close_time, end_to_end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct hdr local_hdr;
    // struct hdr recv_hdr;
    // create_dummy_data(rank, &dummy);
    read_metadata(rank, size, &local_hdr);

    

    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(SOURCE_NAME) - 3;
    // strncpy(filename, SOURCE_NAME, position);
    // strcat(filename, "_new");
    // strcat(filename, SOURCE_NAME + position);
    // if (rank==0) printf("\n%s\n", OUTPUT_NAME);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // MPI_Info_set(info, "nc_hash_size_dim", "4096");
    // MPI_Info_set(info, "nc_hash_size_var", "4096");
    MPI_Info_set(info, "nc_hash_size_dim", "16777216");
    MPI_Info_set(info, "nc_hash_size_var", "8388608"); 
    // MPI_Info_set(info, "nc_hash_size_dim", "2097152");
    // MPI_Info_set(info, "nc_hash_size_var", "1048576"); 
    err = ncmpi_create(MPI_COMM_WORLD, OUTPUT_NAME, cmode, info, &ncid); ERR
    MPI_Barrier(MPI_COMM_WORLD);

    // printf("rank %d, recv_displs: %d, recvcounts: %d \n",  rank, recv_displs[i], recvcounts[i]);

    start_time1 = MPI_Wtime();
    define_hdr(&local_hdr, ncid, rank);
    meta_create_time = MPI_Wtime() - start_time1;
    
    meta_free_hdr(&local_hdr);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time2 = MPI_Wtime();
    err = ncmpi_enddef(ncid); ERR
    enddef_time = MPI_Wtime() - start_time2;

    // Clean up
    MPI_Barrier(MPI_COMM_WORLD);
    // print_memory_info();
    start_time3 = MPI_Wtime();
    err = ncmpi_close(ncid); ERR
    end_time =  MPI_Wtime();
    close_time = end_time - start_time3;
    end_to_end_time = end_time - start_time;
    


    double times[5] = {end_to_end_time, meta_create_time, enddef_time, close_time, total_def_time};
    char *names[5] = {"end-end", "metadata-create", "enddef", "close", "def_dim/var"};
    double max_times[5], min_times[5];

    MPI_Reduce(&times[0], &max_times[0], 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 5; i++) {
        if (rank == 0) {
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
    pnetcdf_check_mem_usage(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;



}