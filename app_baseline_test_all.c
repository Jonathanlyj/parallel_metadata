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

#define FILE_NAME "/files2/scratch/yll6162/parallel_metadata/nue_slice_panoptic_hdf_merged.nc"
// #define FILE_NAME "/files2/scratch/yll6162/parallel_metadata/script/dummy_test.nc"
#define OUTPUT_NAME "app_baseline_test_all.nc"
// #define FILE_NAME "testfile.nc"

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

void read_metdata(int rank, int size, struct hdr *file_info) {

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
            size_t attr_size;
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

/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr(struct hdr *hdr_data, int ncid){
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

static int deserialize_all_hdr(struct hdr **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (struct hdr *)malloc(sizeof(struct hdr));
        deserialize_hdr(all_recv_hdr[i], all_collections_buffer + recv_displs[i], recvcounts[i]);
    }
    return 0;
}

int define_all_hdr(struct hdr **all_recv_hdr, int nproc, int ncid){
    for (int i=0; i< nproc; i++){
        struct hdr *hdr_data = all_recv_hdr[i];
        define_hdr(hdr_data, ncid);
    }
    return 0;
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
    struct hdr dummy;
    // struct hdr recv_hdr;
    // create_dummy_data(rank, &dummy);
    read_metdata(rank, nproc, &dummy);
   // Print the created data for each process
    // printf("\nRank %d:\n", rank);
    // printf("Total Header Size: %lld\n", dummy.xsz);
    // printf("Dimensions:\n");
    // for (int i = 0; i < dummy.dims.ndefined; i++) {
    //     printf("  Name: %s, Size: %lld\n", dummy.dims.value[i]->name, dummy.dims.value[i]->size);
    // }

    // printf("Variables:\n");
    // for (int i = 0; i < dummy.vars.ndefined; i++) {
    //     printf("  Name: %s, Type: %d, NumDims: %d\n", dummy.vars.value[i]->name,  dummy.vars.value[i]->xtype, 
    //     dummy.vars.value[i]->ndims);
    //     printf("    Dim IDs: ");
    //     for (int j = 0; j < dummy.vars.value[i]->ndims; j++) {
    //         printf("%d ", dummy.vars.value[i]->dimids[j]);
    //     }
    //     printf("\n");
    //     printf("    Attributes:\n");
    //     for (int k = 0; k < dummy.vars.value[i]->attrs.ndefined; k++) {
    //         printf("      Name: %s, Nelems: %lld, Type: %d\n", dummy.vars.value[i]->attrs.value[k]->name, 
    //         dummy.vars.value[i]->attrs.value[k]->nelems, dummy.vars.value[i]->attrs.value[k]->xtype);
    //     }
    // }
    // printf("rank %d, buffer size: %lld \n", rank, dummy.xsz);
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = start_time1 = MPI_Wtime();
    char* send_buffer = (char*) malloc(dummy.xsz);
    status = serialize_hdr(&dummy, send_buffer);


    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) malloc(nproc * sizeof(MPI_Offset));
    MPI_Allgather(&dummy.xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

    // Calculate displacements for the second phase
    int* recv_displs = (int*) malloc(nproc * sizeof(int));
    int total_recv_size, min_size, max_size;
    total_recv_size = min_size = max_size = all_collection_sizes[0];
    recv_displs[0] = 0;


    for (int i = 1; i < nproc; ++i) {
        recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
        total_recv_size += all_collection_sizes[i];
        if(all_collection_sizes[i] > max_size){
            max_size = all_collection_sizes[i];
        }
        if(all_collection_sizes[i] < min_size){
            min_size = all_collection_sizes[i];
        }
    }
    double total_recv_size_MB = total_recv_size / (1024.0 * 1024.0);
    double min_size_MB = min_size / (1024.0 * 1024.0);
    double max_size_MB = max_size / (1024.0 * 1024.0);
    // if(rank==0){
    //     printf("\nTotal buffer size: %f MB", total_recv_size_MB);
    //     printf("\nMax buffer size: %f MB", max_size_MB);
    //     printf("\nMin buffer size: %f MB \n", min_size_MB);
    // }
    
    // printf("\nrank %d, dummy xsz %lld", rank, dummy.xsz);
    // Allocate buffer for receiving all header data
    char* all_collections_buffer = (char*) malloc(total_recv_size);
    int* recvcounts =  (int*)malloc(nproc * sizeof(int));
    for (int i = 0; i < nproc; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    // Phase 2: Communicate the actual header data
    // Before MPI_Allgatherv
    MPI_Allgatherv(send_buffer, dummy.xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);
    // Deserialize the received data and print if rank is 0
    
    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(FILE_NAME) - 3;
    // strncpy(filename, FILE_NAME, position);
    // strcat(filename, "_new");
    // strcat(filename, FILE_NAME + position);
    // if (rank==0) printf("\n%s\n", OUTPUT_NAME);
    struct hdr **all_recv_hdr = (struct hdr **)malloc(nproc * sizeof(struct hdr*));
    deserialize_all_hdr(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);

    err = ncmpi_create(MPI_COMM_WORLD, OUTPUT_NAME, cmode, MPI_INFO_NULL, &ncid); ERR
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    define_all_hdr(all_recv_hdr, nproc, ncid);
    // for (int i = 0; i < nproc; ++i) {
    //     struct hdr *recv_hdr = (struct hdr *)malloc(sizeof(struct hdr)); 
    //     deserialize_hdr(recv_hdr, all_collections_buffer + recv_displs[i], recvcounts[i]);
    //     define_hdr(recv_hdr, ncid, rank);
    //     free_hdr(recv_hdr);
    // }
    io_time = MPI_Wtime() - end_time1;
    free_all_hdr(all_recv_hdr, nproc);
    end_time2 = MPI_Wtime();
    err = ncmpi_enddef(ncid); ERR
    end_time3 = MPI_Wtime();
    enddef_time = end_time3 - end_time2;

    // Clean up
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time3 = MPI_Wtime();
    err = ncmpi_close(ncid); ERR
    end_time = MPI_Wtime();
    close_time = end_time - end_time3;
    end_to_end_time = end_time - start_time;
    mpi_time = end_time1 - start_time1;

    double times[6] = {end_to_end_time, mpi_time, io_time, enddef_time, total_def_time, close_time};
    char *names[6] = {"end-end", "mpi-phase", "write", "enddef", "def_dim/var", "close"};
    double max_times[6], min_times[6];


    MPI_Reduce(&times[0], &max_times[0], 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 6, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 6; i++) {
        if (rank == 0) {
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
    pnetcdf_check_mem_usage(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}