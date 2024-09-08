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
#include <mpi.h>
#include <pnetcdf.h>

#define NY 10
#define NX 4

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
    int ncid, varid, blkid, dimid[2];
    char str_att[128];
    float float_att[100];
    MPI_Offset  global_ny, global_nx;
    MPI_Offset start[2], count[2];

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* create a new file for writing ----------------------------------------*/
    cmode |= NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); ERR

    /* the global array is NY * (NX * nprocs) */
    global_ny = NY;
    global_nx = NX * (rank + 1);


    /* add a global attribute: a time stamp at rank 0 */
    time_t ltime = time(NULL); /* get the current calendar time */
    asctime_r(localtime(&ltime), str_att);
    sprintf(str_att, "Mon Aug 13 21:27:48 2018");


    // if (rank < 3){
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "history", strlen(str_att),
                             &str_att[0]); ERR

    /* define dimensions x and y */
    char str_y[20];
    char str_blk[20];
    // char shared_y[20];
    // char shared_x[20];
    char str_x[20];
    char var_rank[20];
    char shared_var[20];
    // sprintf(shared_y, "Y_shared");
    // sprintf(shared_x, "X_shared");
    sprintf(str_y, "dim_Y_rank_%d", rank);
    sprintf(str_x, "dim_X_rank_%d", rank);
    sprintf(str_blk, "blk_rank_%d", rank);

    err = ncmpi_def_block(ncid, str_blk, &blkid); ERR
    err = ncmpi_def_dim(ncid, blkid, str_y, NY, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, blkid, str_x, global_nx, &dimid[1]); ERR
    int buf[NY][global_nx];
    for (i=0; i<NY; i++)
        for (j=0; j<global_nx; j++)
            buf[i][j] = rank + 1;
    sprintf(var_rank, "var_rank_%d", rank);

    err = ncmpi_def_var(ncid, blkid, var_rank, NC_INT, 2, dimid, &varid); ERR
    err = ncmpi_enddef(ncid); ERR


    err = ncmpi_put_var_int_all(ncid, blkid, varid,  &buf[0][0]); ERR
    


    MPI_Barrier(comm);
    err = ncmpi_close(ncid); ERR

    // printf("\n rank %d: after file close", rank);
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
    if (argv[optind] == NULL) strcpy(filename, "new_format_create_simple.pnc");
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