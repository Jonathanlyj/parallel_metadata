/*********************************************************************
 *
 *  Copyright (C) 2023, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/
/* $Id$ */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    To compile:
 *       mpicc -O2 no_comm.c -o no_comm -lpnetcdf
 *
 * Example commands for MPI run and outputs from running ncmpidump on the
 * output netCDF file produced by this example program:
 * 
 *  mpiexec -n 4 ./no_comm

netcdf no_comm {
// file format: CDF-1
dimensions:
        dim_0_1 = 1 ;
        dim_0_0 = 1 ;
        dim_1_1 = 2 ;
        dim_1_0 = 2 ;
        dim_2_1 = 3 ;
        dim_2_0 = 3 ;
        dim_3_1 = 4 ;
        dim_3_0 = 4 ;
variables:
        int var_0_0(dim_0_1, dim_0_0) ;
        int var_1_0(dim_1_1, dim_1_0) ;
        int var_2_0(dim_2_1, dim_2_0) ;
        int var_3_0(dim_3_1, dim_3_0) ;
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
    
    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); ERR

    int *varid = (int *)malloc(nprocs * sizeof(int));
    int *dimid = (int *)malloc(2 * nprocs * sizeof(int));



    for (i=0; i<nprocs; i++){
        /* define dimensions x and y */
        sprintf(dim_x_name, "dim_%d_0", i);
        sprintf(dim_y_name, "dim_%d_1", i);
        sprintf(var_name, "var_%d_0", i);
        // sprintf(float_att_name, "float_att_name_%d", i);
        err = ncmpi_def_dim(ncid, dim_y_name, i + 1, &dimid[2 * i]); ERR
        err = ncmpi_def_dim(ncid, dim_x_name, i + 1, &dimid[2 * i + 1]); ERR

        /* define a 2D variable of integer type */
        err = ncmpi_def_var(ncid, var_name, NC_INT, 2,  &dimid[2 * i], &varid[i]); ERR

        // for (j=0; j<8; j++) float_att[j] = j;
        // err = ncmpi_put_att_float(ncid, varid[i], float_att_name, NC_FLOAT, 8,
        //                         &float_att[0]); ERR
    }


    int buf[rank + 1][rank + 1];
    for (i=0; i<rank + 1; i++)
        for (j=0; j<rank + 1; j++)
             buf[i][j] = rank;
    /* do not forget to exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* now we are in data mode */

    err = ncmpi_put_var_int_all(ncid, varid[rank], &buf[0][0]); ERR

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

    if (verbose && rank == 0) printf("%s: example of using put_vara APIs\n",__FILE__);

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