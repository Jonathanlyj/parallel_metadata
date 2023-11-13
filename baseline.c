/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    To compile:
 *       mpicc -O2 baseline.c -o baseline -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output netCDF file produced by this example program:
 * 
 *  mpiexec -n 4 ./baseline ./baseline.nc

netcdf baseline {
// file format: CDF-1
dimensions:
        dim_0_0 = 1 ;
        dim_0_1 = 1 ;
        dim_1_0 = 2 ;
        dim_1_1 = 2 ;
        dim_2_0 = 3 ;
        dim_2_1 = 3 ;
        dim_3_0 = 4 ;
        dim_3_1 = 4 ;
variables:
        int var_0_0(dim_0_0, dim_0_1) ;
        int var_1_0(dim_1_0, dim_1_1) ;
        int var_2_0(dim_2_0, dim_2_1) ;
        int var_3_0(dim_3_0, dim_3_1) ;
}
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <mpi.h>
#include <pnetcdf.h>



static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] | [-q] [-k format] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-k format] file format: 1 for CDF-1, 2 for CDF-2, 3 for NetCDF4,\n"
    "                                4 for NetCDF4 classic model, 5 for CDF-5\n"
    "       [filename] output netCDF file name\n";
    fprintf(stderr, help, argv0);
}

/*----< pnetcdf_io() >-------------------------------------------------------*/
static int
pnetcdf_io(MPI_Comm comm, char *filename, int cmode)
{
    int i, j, rank, nprocs, err, nerrs=0;
    int ncid;
    char str_att[128];
    float float_att[100];
    char dim_x_name[10];
    char dim_y_name[10];
    char var_name[10];
    char float_att_name[10];

    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int buf[rank + 1][rank + 1];

    // int *varid = (int *)malloc(nprocs * sizeof(int));
    // int *dimid = (int *)malloc(2 * nprocs * sizeof(int));

    // Each process has a different number of strings, but all strings have the same length
    int total_dim_count = 0;
    int local_dim_count = 2;
    int dim_name_len = 10; 

    int total_var_count = 0, total_var_dim_num = 0;
    int local_var_count = 1;
    int var_name_len = 20; 




    /* ----------------------------------  DIMENSIONs  ----------------------------------------*/
    char local_dims_names[local_dim_count][dim_name_len];
    int local_dims_sizes[local_dim_count];
    // Initialize local dimension names and sizes 
    for (int i = 0; i < local_dim_count; ++i) {
        sprintf(local_dims_names[i], "dim_%d_%d", rank, i);
        local_dims_sizes[i] = rank + 1;
    }



    // Gather the number of dims to be sent to each process
    int dim_counts[nprocs];
    int dim_char_counts[nprocs];
    int displs_dim[nprocs];

    MPI_Allgather(&local_dim_count, 1, MPI_INT, dim_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    for(int i=0; i<nprocs; i++) total_dim_count += dim_counts[i];
    // Exchange dim names data
    // Calculate displacement values for MPI_Alltoallv
    displs_dim[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        displs_dim[i] = displs_dim[i - 1] + dim_counts[i - 1] * dim_name_len;
    }
    // Allocate memory for receiving dim names from all processes
    char all_dim_names[total_dim_count][dim_name_len]; 

    for(int i=0; i<nprocs; i++) dim_char_counts[i] = dim_counts[i] * dim_name_len;

    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_dims_names, dim_char_counts[rank], MPI_CHAR, all_dim_names,
                  dim_char_counts, displs_dim, MPI_CHAR, MPI_COMM_WORLD);

    // Exchange dim size data
    // Calculate displacement values for MPI_Alltoallv
    displs_dim[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        displs_dim[i] = displs_dim[i - 1] + dim_counts[i - 1];
    }
    // Allocate memory for receiving dim names from all processes
    int all_dim_sizes[total_dim_count]; 
    for(int i=0; i<nprocs; i++) dim_char_counts[i] = dim_counts[i] * dim_name_len;
    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_dims_sizes, dim_counts[rank], MPI_INT, all_dim_sizes,
                  dim_counts, displs_dim, MPI_INT, MPI_COMM_WORLD);
    // // Print the received dims
    // for (int i = 0; i < total_dim_count; ++i) {
    //     printf("Rank %d received: %s size: %d\n", rank, all_dim_names[i], all_dim_sizes[i]);
        
    // }

    /* ----------------------------------  VARIABLES  ----------------------------------------*/
    char local_vars_names[local_var_count][var_name_len];
    int local_var_num_dim[local_var_count];
    int local_var_type[local_var_count];
    int local_var_dim_idx[2 * local_var_count];
    // Initialize local variaable names, num_dims, var type, and dim_ids
    for (int i = 0; i < local_var_count; ++i) {
        sprintf(local_vars_names[i], "var_%d_%d", rank, i);
        local_var_num_dim[i] = 2;
        local_var_type[i] = NC_INT;
        for (int j=0; j < local_var_num_dim[i]; ++j)
            local_var_dim_idx[2 * i + j] = j;
    }
    // Gather the number of vars to be sent to each process
    int var_counts[nprocs];
    int var_char_counts[nprocs];
    int displs_var[nprocs];

    MPI_Allgather(&local_var_count, 1, MPI_INT, var_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    for(int i=0; i<nprocs; i++) total_var_count += var_counts[i];
    // Exchange var names data
    // Calculate displacement values for MPI_Alltoallv
    displs_var[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        displs_var[i] = displs_var[i - 1] + var_counts[i - 1] * var_name_len;
    }
    // Allocate memory for receiving var names from all processes
    char all_var_names[total_var_count][var_name_len]; 

    for(int i=0; i<nprocs; i++) var_char_counts[i] = var_counts[i] * var_name_len;

    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_vars_names, var_char_counts[rank], MPI_CHAR, all_var_names,
                  var_char_counts, displs_var, MPI_CHAR, MPI_COMM_WORLD);

    // Exchange var num_dim data
    // Calculate displacement values for MPI_Alltoallv
    displs_var[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        displs_var[i] = displs_var[i - 1] + var_counts[i - 1];
    }
    // Allocate memory for receiving var num_dims from all processes
    int all_var_num_dims[total_var_count]; 
    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_var_num_dim, var_counts[rank], MPI_INT, all_var_num_dims,
                  var_counts, displs_var, MPI_INT, MPI_COMM_WORLD);
    for(int i=0; i < total_var_count; ++i) total_var_dim_num += all_var_num_dims[i];
    // Exchange var type data
    // Allocate memory for receiving var names from all processes
    int all_var_type[total_var_count]; 
    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_var_type, var_counts[rank], MPI_INT, all_var_type,
                  var_counts, displs_var, MPI_INT, MPI_COMM_WORLD);

    // Exchange var dim index data
    int displs_var_dim_idx[nprocs];
    int all_var_dim_idx[total_var_dim_num];
    int var_dim_idx_count[nprocs];
    // Calculate count and displacement values for MPI_Alltoallv
    
    for (int i = 0; i < nprocs; ++i) {
        var_dim_idx_count[i] = 0;
        for (int j = 0; j < var_counts[i]; ++j)
            var_dim_idx_count[i] += all_var_num_dims[displs_var[i] + j];
    }
    displs_var_dim_idx[0] = 0;
    for (int i = 1; i < nprocs; ++i) 
        displs_var_dim_idx[i] = displs_var_dim_idx[i - 1] + var_dim_idx_count[i - 1];
    
    // Use MPI_Alltoallv to exchange data
    MPI_Allgatherv(local_var_dim_idx, 2 * local_var_count, MPI_INT, all_var_dim_idx,
                  var_dim_idx_count, displs_var_dim_idx, MPI_INT, MPI_COMM_WORLD);
                
    // Print the received vars
    // int num_dim_idx_ct;
    // for (int i = 0; i < nprocs; ++i){
    //     num_dim_idx_ct = 0;
    //     for (int j = 0; j < var_counts[i]; j++){
    //         printf("Rank %d received: %s num_dim: %d, var_type: %d\n", rank, all_var_names[displs_var[i] + j], all_var_num_dims[displs_var[i] + j], all_var_type[displs_var[i] + j]);
    //         for (int k = 0; k < all_var_num_dims[displs_var[i] + j]; k++)
    //             printf("     dim %d name: %s\n", k, all_dim_names[displs_dim[i] + all_var_dim_idx[displs_var_dim_idx[i] + num_dim_idx_ct + k]]);
    //         num_dim_idx_ct += all_var_num_dims[displs_var[i] + j];
    //     }
    // }



    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); ERR

    int *varid = (int *)malloc(total_var_count * sizeof(int));
    int *dimid = (int *)malloc(total_dim_count * sizeof(int));

    for (i=0; i<rank + 1; i++)
        for (j=0; j<rank + 1; j++)
             buf[i][j] = rank;
    
    int dim_idx, var_idx, num_dim_idx_ct, var_dim_id_idx;
    int local_var_ids[local_var_count];
    for (i=0; i<nprocs; i++){
        /* define local dimensions for this process */
        // sprintf(dim_x_name, "rank_%d_x", i);
        // sprintf(dim_y_name, "rank_%d_y", i);
        // sprintf(var_name, "var_%d", i);
        // sprintf(float_att_name, "float_att_name_%d", i);
        for (j=0; j < dim_counts[i]; ++j){
            dim_idx = displs_dim[i] + j;
            err = ncmpi_def_dim(ncid, all_dim_names[dim_idx], all_dim_sizes[dim_idx], &dimid[dim_idx]); ERR
        }
        

        /* define a 2D variable of integer type */
        num_dim_idx_ct = 0;
        for (j=0; j < var_counts[i]; ++j){
            var_idx = displs_var[i] + j;
            var_dim_id_idx = displs_dim[i] + all_var_dim_idx[displs_var_dim_idx[i] + num_dim_idx_ct];
            err = ncmpi_def_var(ncid, all_var_names[var_idx], all_var_type[var_idx], all_var_num_dims[var_idx],  &dimid[var_dim_id_idx], &varid[var_idx]); ERR
            if (i == rank) local_var_ids[j] = varid[var_idx];
            num_dim_idx_ct += all_var_num_dims[displs_var[i] + j];
        }
        

    }



    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* now we are in data mode */
    for (i = 0; i < local_var_count; i++){
        err = ncmpi_put_var_int_all(ncid, local_var_ids[i], &buf[0][0]); ERR
    }


    err = ncmpi_close(ncid); ERR

    return nerrs;
}

int main(int argc, char** argv)
{
    extern int optind;
    extern char *optarg;
    char filename[256];
    int i, rank, kind=0, cmode=0, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqk:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': kind = atoi(optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }
    if (argv[optind] == NULL) strcpy(filename, "testfile.nc");
    else                      snprintf(filename, 256, "%s", argv[optind]);

    MPI_Bcast(filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD);


    switch (kind) {
        case(2): cmode = NC_64BIT_OFFSET;             break;
        case(3): cmode = NC_NETCDF4;                  break;
        case(4): cmode = NC_NETCDF4|NC_CLASSIC_MODEL; break;
        case(5): cmode = NC_64BIT_DATA;               break;
        default: cmode = 0;
    }

#ifndef PNETCDF_DRIVER_NETCDF4
    /* netcdf4 driver is not enabled, skip */
    if (kind == 3 || kind == 4) {
        MPI_Finalize();
        return 0;
    }
#endif
    clock_t start_time = clock();
    nerrs += pnetcdf_io(MPI_COMM_WORLD, filename, cmode);
    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time) / CLOCKS_PER_SEC);
    printf("CPU time used: %f seconds\n", cpu_time_used);
    MPI_Finalize();
    return (nerrs > 0);
}