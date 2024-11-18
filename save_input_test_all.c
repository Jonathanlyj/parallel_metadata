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
#include "baseline_ncx_app.h" 



static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

// #define FILE_NAME "/global/homes/y/yll6162/parallel_metadata/data/nue_slice_panoptic_hdf_merged.nc"
// #define FILE_NAME "/files2/scratch/yll6162/parallel_metadata/script/nue_slice_panoptic_hdf_merged_10_copy.nc"
#define FILE_NAME "/global/homes/y/yll6162/parallel_metadata/data/nue_slice_panoptic_hdf_merged_10_copy.nc"
// #define FILE_NAME "/files2/scratch/yll6162/parallel_metadata/script/local_hdr_test.nc"
// #define OUTPUT_NAME "/pscratch/sd/y/yll6162/FS_2M_32/save_input_test_all"
// #define FILE_NAME "testfile.nc"
#define OUTPUT_NAME "/pscratch/sd/y/yll6162/FS_2M_8/save_input_test_all_10_copy"


double def_start_time, total_def_time=0;
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


/* ---------------------------------- Read Metadata ----------------------------------------*/

void read_metadata(int rank, int nproc, struct hdr *file_info) {

    int ncid, num_vars, num_dims, tot_num_dims, elem_sz, v_attrV_xsz, status;
    MPI_Offset start, count;
    file_info->dims.ndefined = 0;
    file_info->attrs.ndefined = 0;
    file_info->vars.ndefined = 0;
    file_info->xsz = 0;
    tot_num_dims = 0;
    // Open the NetCDF file
    if (ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, MPI_INFO_NULL, &ncid) != NC_NOERR) {
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
    int vars_per_process = num_vars / nproc;
    int remainder = num_vars % nproc;
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
            xlen_nc_type(attr_type, &elem_sz);
            variable_info->attrs.value[j]->xvalue = malloc(attr_size * elem_sz);
            ncmpi_get_att(ncid, i, variable_info->attrs.value[j]->name, variable_info->attrs.value[j]->xvalue);

            file_info->xsz += sizeof(uint32_t) + sizeof(char) * variable_info->attrs.value[j]->name_len; //attr name
            file_info->xsz += sizeof(uint32_t); // nc_type
            file_info->xsz += sizeof(uint32_t); // nelems
            status = xlen_nc_type(variable_info->attrs.value[j]->xtype, &v_attrV_xsz);
            file_info->xsz += variable_info->attrs.value[j]->nelems * v_attrV_xsz; // attr_value
        }

        // Add the variable information to the file_info structure
        file_info->vars.value[i - start] = variable_info;
    }

    // Close the NetCDF file
    ncmpi_close(ncid);
    pnetcdf_check_mem_usage(MPI_COMM_WORLD);
}




static int free_all_hdr(struct hdr **all_recv_hdr, int nproc){
    if (all_recv_hdr != NULL){
        for (int i=0; i< nproc; i++) free_hdr(all_recv_hdr[i]);
        free(all_recv_hdr);
    }
    return 0;
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, nerrs=0;
    double end_to_end_time, mpi_time, io_time, enddef_time, close_time, max_time, min_time;
    double start_time, start_time1, end_time1, end_time2, end_time3, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    struct hdr local_hdr;
    // struct hdr recv_hdr;
    // create_local_hdr_data(rank, &local_hdr);
    read_metadata(rank, nproc, &local_hdr);
   // Print the created data for each process
    // printf("\nRank %d:\n", rank);
    // printf("Total Header Size: %lld\n", local_hdr.xsz);
    // printf("Dimensions:\n");
    // for (int i = 0; i < local_hdr.dims.ndefined; i++) {
    //     printf("  Name: %s, Size: %lld\n", local_hdr.dims.value[i]->name, local_hdr.dims.value[i]->size);
    // }

    // printf("Variables:\n");
    // for (int i = 0; i < local_hdr.vars.ndefined; i++) {
    //     printf("  Name: %s, Type: %d, NumDims: %d\n", local_hdr.vars.value[i]->name,  local_hdr.vars.value[i]->xtype, 
    //     local_hdr.vars.value[i]->ndims);
    //     printf("    Dim IDs: ");
    //     for (int j = 0; j < local_hdr.vars.value[i]->ndims; j++) {
    //         printf("%d ", local_hdr.vars.value[i]->dimids[j]);
    //     }
    //     printf("\n");
    //     printf("    Attributes:\n");
    //     for (int k = 0; k < local_hdr.vars.value[i]->attrs.ndefined; k++) {
    //         printf("      Name: %s, Nelems: %lld, Type: %d\n", local_hdr.vars.value[i]->attrs.value[k]->name, 
    //         local_hdr.vars.value[i]->attrs.value[k]->nelems, local_hdr.vars.value[i]->attrs.value[k]->xtype);
    //     }
    // }
    // printf("rank %d, buffer size: %lld \n", rank, local_hdr.xsz);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = start_time1 = MPI_Wtime();
    char* send_buffer = (char*) malloc(local_hdr.xsz);
    printf("\nlocal_hdr.xsz: %lld\n", local_hdr.xsz);
    status = serialize_hdr(&local_hdr, send_buffer);

    FILE* file = fopen(OUTPUT_NAME, "wb");

    if (file == NULL) {
        perror("Failed to open file");
        goto fn_exit;
        return 1;
    }
    // Write the buffer to the file
    size_t written = fwrite(send_buffer, 1, local_hdr.xsz, file);
    if (written != local_hdr.xsz) {
        perror("Failed to write complete buffer to file");
    }
    fclose(file);

    // for (int i = 0; i < nproc; ++i) {
    //     struct hdr *recv_hdr = (struct hdr *)malloc(sizeof(struct hdr)); 
    //     deserialize_hdr(recv_hdr, all_collections_buffer + recv_displs[i], recvcounts[i]);
    //     define_hdr(recv_hdr, ncid, rank);
    //     free_hdr(recv_hdr);
    // }


    // Clean up
fn_exit:
    free(send_buffer);

    MPI_Barrier(MPI_COMM_WORLD);



    MPI_Finalize();
    return 0;
}