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
typedef struct {
    int total_names_length;
    int num_dimensions;
    int* dimension_sizes;
    int* name_lengths;
    char* dimension_names;
} DimensionCollection;

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



// Function to serialize DimensionCollection into a byte buffer
void serialize(DimensionCollection* collection, char* buffer) {
    int offset = 0;

    // Copy total_names_length and num_dimensions
    memcpy(buffer + offset, &(collection->total_names_length), sizeof(int));
    offset += sizeof(int);
    memcpy(buffer + offset, &(collection->num_dimensions), sizeof(int));
    offset += sizeof(int);

    // Copy dimension_sizes and name_lengths arrays
    memcpy(buffer + offset, collection->dimension_sizes, collection->num_dimensions * sizeof(int));
    offset += collection->num_dimensions * sizeof(int);
    memcpy(buffer + offset, collection->name_lengths, collection->num_dimensions * sizeof(int));
    offset += collection->num_dimensions * sizeof(int);

    // Copy dimension_names string
    memcpy(buffer + offset, collection->dimension_names, collection->total_names_length * sizeof(char));
}

// Function to deserialize a byte buffer into DimensionCollection
void deserialize(char* buffer, DimensionCollection* collection) {
    int offset = 0;

    // Copy total_names_length and num_dimensions
    memcpy(&(collection->total_names_length), buffer + offset, sizeof(int));
    offset += sizeof(int);
    memcpy(&(collection->num_dimensions), buffer + offset, sizeof(int));
    offset += sizeof(int);

    // Allocate memory for dimension_sizes and name_lengths arrays
    collection->dimension_sizes = (int*) malloc(collection->num_dimensions * sizeof(int));
    collection->name_lengths = (int*) malloc(collection->num_dimensions * sizeof(int));

    // Copy dimension_sizes and name_lengths arrays
    memcpy(collection->dimension_sizes, buffer + offset, collection->num_dimensions * sizeof(int));
    offset += collection->num_dimensions * sizeof(int);
    memcpy(collection->name_lengths, buffer + offset, collection->num_dimensions * sizeof(int));
    offset += collection->num_dimensions * sizeof(int);

    // Allocate memory and copy dimension_names string
    collection->dimension_names = (char*) malloc(collection->total_names_length * sizeof(char));
    memcpy(collection->dimension_names, buffer + offset, collection->total_names_length * sizeof(char));
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


    // Example: Each process has a number of dimensions equal to its rank + 1
    int num_dimensions = rank + 1;

    // Allocate and fill dimension collection for the current process
    DimensionCollection send_collection;
    send_collection.num_dimensions = num_dimensions;
    send_collection.dimension_sizes = (int*) malloc(num_dimensions * sizeof(int));
    send_collection.name_lengths = (int*) malloc(num_dimensions * sizeof(int));
    send_collection.total_names_length = 0;
    char* names_concat = (char*) malloc(num_dimensions * 50 * sizeof(char)); // Each name up to 50 chars
    names_concat[0] = '\0';
    for (int i = 0; i < num_dimensions; ++i) {
        char dimension_name[50];
        sprintf(dimension_name, "Process_%d_Dim_%d", rank, i);
        send_collection.dimension_sizes[i] = rank * 10 + i;
        send_collection.name_lengths[i] = strlen(dimension_name) + 1;
        send_collection.total_names_length += send_collection.name_lengths[i];
        strcat(names_concat, dimension_name);
        strcat(names_concat, ";"); // Separator
    }
    send_collection.dimension_names = names_concat;

    // Serialize the send_collection
    int my_collection_size = sizeof(int) * (2 + 2 * send_collection.num_dimensions) + sizeof(char) * send_collection.total_names_length;
    char* send_buffer = (char*) malloc(my_collection_size);
    serialize(&send_collection, send_buffer);

    // Phase 1: Communicate the sizes of the DimensionCollection structure for each process
    int* all_collection_sizes = (int*) malloc(nprocs * sizeof(int));
    MPI_Allgather(&my_collection_size, 1, MPI_INT, all_collection_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    // Calculate displacements for the second phase
    int* recv_displs = (int*) malloc(nprocs * sizeof(int));
    int total_recv_size = all_collection_sizes[0];
    recv_displs[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
        recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
        total_recv_size += all_collection_sizes[i];
    }

    // Allocate buffer for receiving all DimensionCollection data
    char* all_collections_buffer = (char*) malloc(total_recv_size);

    // Phase 2: Communicate the actual DimensionCollection data
    MPI_Allgatherv(send_buffer, my_collection_size, MPI_BYTE,
                   all_collections_buffer, all_collection_sizes, recv_displs, MPI_BYTE, MPI_COMM_WORLD);

    // Deserialize the received data
    // For example, deserialize the data for the first process (process 0)

// Deserialize the received data and print if rank is 0
if (rank == 0){
    for (int i = 0; i < nprocs; ++i) {
        DimensionCollection received_collection;
        deserialize(all_collections_buffer + recv_displs[i], &received_collection);

        printf("Process %d:\n", i);
        printf("  Total Name Length: %d\n", received_collection.total_names_length);
        printf("  Number of Dimensions: %d\n", received_collection.num_dimensions);
        char* name_ptr = received_collection.dimension_names;
        for (int j = 0; j < received_collection.num_dimensions; ++j) {
            printf("    Dimension %s: Size = %d, Name Length = %d\n",
                   name_ptr, received_collection.dimension_sizes[j], received_collection.name_lengths[j]);
            name_ptr += received_collection.name_lengths[j]; 
        }
        
        // Free the arrays inside received_collection
        free(received_collection.dimension_sizes);
        free(received_collection.name_lengths);
        free(received_collection.dimension_names);
    }
}
    // ... Process the deserialized data ...

    // Clean up
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    free(send_collection.dimension_sizes);
    free(send_collection.name_lengths);
    free(send_collection.dimension_names);
    // ... Also, free the arrays inside received_collection ...

    // MPI_Finalize();
    return 0;




    // /* create a new file for writing ----------------------------------------*/
    // cmode |= NC_CLOBBER;
    // err = ncmpi_create(comm, filename, cmode, MPI_INFO_NULL, &ncid); ERR

    // int *varid = (int *)malloc(total_var_count * sizeof(int));
    // int *dimid = (int *)malloc(total_dim_count * sizeof(int));

    // for (i=0; i<rank + 1; i++)
    //     for (j=0; j<rank + 1; j++)
    //          buf[i][j] = rank;
    
    // int dim_idx, var_idx, num_dim_idx_ct, var_dim_id_idx;
    // int local_var_ids[local_var_count];
    // for (i=0; i<nprocs; i++){
    //     /* define local dimensions for this process */
    //     // sprintf(dim_x_name, "rank_%d_x", i);
    //     // sprintf(dim_y_name, "rank_%d_y", i);
    //     // sprintf(var_name, "var_%d", i);
    //     // sprintf(float_att_name, "float_att_name_%d", i);
    //     for (j=0; j < dim_counts[i]; ++j){
    //         dim_idx = displs_dim[i] + j;
    //         err = ncmpi_def_dim(ncid, all_dim_names[dim_idx], all_dim_sizes[dim_idx], &dimid[dim_idx]); ERR
    //     }
        

    //     /* define a 2D variable of integer type */
    //     num_dim_idx_ct = 0;
    //     for (j=0; j < var_counts[i]; ++j){
    //         var_idx = displs_var[i] + j;
    //         var_dim_id_idx = displs_dim[i] + all_var_dim_idx[displs_var_dim_idx[i] + num_dim_idx_ct];
    //         err = ncmpi_def_var(ncid, all_var_names[var_idx], all_var_type[var_idx], all_var_num_dims[var_idx],  &dimid[var_dim_id_idx], &varid[var_idx]); ERR
    //         if (i == rank) local_var_ids[j] = varid[var_idx];
    //         num_dim_idx_ct += all_var_num_dims[displs_var[i] + j];
    //     }
        

    // }



    // /* do not forget to exit define mode */
    // err = ncmpi_enddef(ncid); ERR

    // /* now we are in data mode */
    // for (i = 0; i < local_var_count; i++){
    //     err = ncmpi_put_var_int_all(ncid, local_var_ids[i], &buf[0][0]); ERR
    // }


    // err = ncmpi_close(ncid); ERR

    // return nerrs;
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