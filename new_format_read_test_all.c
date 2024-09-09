/*********************************************************************
 *
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <time.h>   /* time() localtime(), asctime() */
#include <assert.h>
#include <mpi.h>
#include <pnetcdf.h>


#define NY 10
#define NX 4

static int verbose;

#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}
#define SOURCE_NAME "/pscratch/sd/y/yll6162/FS_2M_32/new_format_test_all_256.pnc"
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
    int ncid, varid, blkid, dimid[2];
    char str_att[128];
    float float_att[100];

    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    
    double total_read_time = 0;
    double read_start_time = MPI_Wtime();
    /* open the newly created file for read only -----------------------------*/
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    ERR
    /*read the number of blocks*/
    int nblocks;
    err = ncmpi_inq_nblocks(ncid, &nblocks);
    // printf("\nThere are this number of blocks: %d\n", nblocks);

    /* read the block name & block id*/
    char blk_name[20], var_name[20];    
    int nvars;
    int ndims, v_ndims;
    int* v_dimids;
    MPI_Offset* dim_sizes;
    MPI_Offset var_size = 1;


    for (int blkid = 0; blkid < nblocks; blkid++){
        // err = ncmpi_inq_blockname(ncid, blkid, blk_name);
        // printf("\nBlock %d name: %s\n", blkid, blk_name);
        // err = ncmpi_inq_blkid(ncid, blk_name, &blkid);
        err= ncmpi_open_block(ncid, blkid);
        ERR;
        err = ncmpi_inq_block(ncid, blkid, NULL, &ndims, &nvars, NULL);
        ERR;
        
        for (varid = 0; varid < nvars; varid++){
            err = ncmpi_inq_var(ncid, blkid, varid, var_name, NULL, &v_ndims, NULL, NULL);
            // printf("\nBlock %d has variable %s with %d dims\n", blkid, var_name, v_ndims);
            v_dimids = (int*) malloc(v_ndims * sizeof(int));
            dim_sizes = (MPI_Offset*) malloc(v_ndims * sizeof(MPI_Offset));
            err = ncmpi_inq_var(ncid, blkid, varid, NULL, NULL, NULL, v_dimids, NULL);
            for (int j = 0; j < v_ndims; j++){
                err = ncmpi_inq_dimlen(ncid, blkid, v_dimids[j], &dim_sizes[j]);
                // var_size *= dim_sizes[j];
                // printf("\nDimension %d has size %lld\n", j, dim_sizes[j]);
            }
            free(dim_sizes);
            free(v_dimids);
        }


    }
    // total_read_time += MPI_Wtime() - read_start_time;


    

    // printf("\nRank %d will read block %d\n", rank, blkid);

    // printf("\nBlock %d has %d dimensions and %d variables\n", blkid, ndims, nvars);

    // Allocate memory for the raw data array (since we're dealing with integers)
    // printf("\nVariable size: %lld\n", var_size);
    // int *data = (int *)malloc(var_size * sizeof(int));
    // if (data == NULL) {
    //     printf("Error: memory allocation failed\n");
    //     MPI_Abort(MPI_COMM_WORLD, 1);
    // }
    // // Read the integer variable data into the C array
    // if ((err = ncmpi_get_var_int_all(ncid, blkid, 0, data))) {
    //     printf("Error: ncmpi_get_var_int_all() failed (%s)\n", ncmpi_strerror(err));
    //     MPI_Abort(MPI_COMM_WORLD, err);
    // }

    // Print a small portion of the data (for example, the first 10 values)
    // printf("First 10 data values: ");
    // for (int i = 0; i < 10 && i < var_size; i++) {
    //     printf("%d ", data[i]);
    // }
    // ERR;
    // free(data);


    /* close file */
    err = ncmpi_close(ncid);
    ERR
    total_read_time += MPI_Wtime() - read_start_time;
    double max_time, min_time;
    MPI_Reduce(&total_read_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_read_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Max read time: %f seconds\n",  max_time);
        printf("Min read time: %f seconds\n", min_time);
    }

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
    if (argv[optind] == NULL) strcpy(filename, SOURCE_NAME);
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
    nerrs += pnetcdf_io(MPI_COMM_WORLD, filename, cmode);

    MPI_Finalize();
    return (nerrs > 0);
}