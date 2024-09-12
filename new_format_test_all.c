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
#include <assert.h>
// #include "baseline_ncx_lib.h"
#include "baseline_ncx_app.h" 




static int verbose;



#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}

// #define SOURCE_NAME "/pscratch/sd/y/yll6162/FS_2M_32/save_input_test_all"
#define SOURCE_NAME "save_input_test_all_10_copy"
// #define SOURCE_NAME "/files2/scratch/yll6162/parallel_metadata/script/dummy_test.nc"
// #define OUTPUT_NAME "/pscratch/sd/y/yll6162/FS_2M_32/new_format_test_all.pnc"
#define OUTPUT_NAME "new_format_test_all_10_copy.pnc"

double def_start_time;
double total_def_time = 0;


//Following functions definitions are going to be overridden by LD_PRELOAD
// int ncmpi_create(MPI_Comm comm, const char *path, int cmode, MPI_Info info, int *ncidp){
//     printf("this function is not supposed to run\n");
//     return 0;
// }

// int ncmpi_enddef(int ncid) {
//     printf("this function is not supposed to run\n");
//     return 0;
// }

// int ncmpi_close(int ncid){
//     printf("this function is not supposed to run\n");
//     return 0;
// }

// int ncmpi_def_block(int ncid, const char *name, int *blkidp){
//     printf("this function is not supposed to run\n");
//     return 0;
// }

// int ncmpi_def_var(int ncid, int blkid, const char *name, nc_type xtype, int ndims, const int *dimidsp, int *varidp){
//     printf("this function is not supposed to run\n");
//     return 0;
// }

// int ncmpi_def_dim(int ncid, int blkid, const char *name, MPI_Offset len, int *idp){
//     printf("this function is not supposed to run\n");
//     return 0;
// }

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

int read_metadata_from_file(const char* filename, struct hdr *recv_hdr) {
    // Open the file in binary read mode
    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

    // Move file pointer to the end to determine the file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);  // Move the file pointer back to the beginning

    // Allocate memory for the buffer
    char *buffer = (char*) malloc(file_size);
    if (buffer == NULL) {
        perror("Failed to allocate memory");
        fclose(file);
        return 1;
    }

    // Read the file contents into the buffer
    // printf("file_size: %ld\n", file_size);
    size_t read_size = fread(buffer, 1, file_size, file);
    if (read_size != file_size) {
        perror("Failed to read complete file");
        free(buffer);
        fclose(file);
        return 1;
    }

    // Close the file after reading
    fclose(file);
    recv_hdr->xsz = file_size;
    deserialize_hdr(recv_hdr, buffer, file_size);
    free(buffer);
    return 0;
}
  

/* ---------------------------------- Decode Metadata ----------------------------------------*/
int define_hdr_nf(struct hdr *hdr_data, int ncid, int rank, int var_start, int var_count) {
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    
    int i,j,k,nerrs=0;
    int err;
    def_start_time = MPI_Wtime();
    char str_blk[20];
    sprintf(str_blk, "blk_rank_%d", rank);
    int blkid;

    def_start_time = MPI_Wtime();
    err = ncmpi_def_block(ncid, str_blk, &blkid);
    total_def_time += MPI_Wtime() - def_start_time;

    //define variables
    int *varid = (int *)malloc(var_count * sizeof(int));
    int v_ndims, v_namelen, xtype, n_att;
    int att_namelen, att_xtype, att_nelems;

    for (i=var_start; i<var_start+var_count; i++){
        
        int varindex = i - var_start;
        v_namelen =  hdr_data->vars.value[i]->name_len;
        xtype = hdr_data->vars.value[i]->xtype;

        v_ndims = hdr_data->vars.value[i]->ndims;
        int *v_dimids = (int *)malloc(v_ndims * sizeof(int));
        for (j=0; j<v_ndims; j++){
            int dimindex = hdr_data->vars.value[i]->dimids[j];
            def_start_time = MPI_Wtime();
            err = ncmpi_def_dim(ncid, blkid, hdr_data->dims.value[dimindex]->name,  hdr_data->dims.value[dimindex]->size, &v_dimids[j]); ERR
            total_def_time += MPI_Wtime() - def_start_time;
            // }
        }

        def_start_time = MPI_Wtime();

        if (v_ndims == 0)
            err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, NULL, &varid[varindex]);  
        else
            err = ncmpi_def_var(ncid, blkid, hdr_data->vars.value[i]->name, xtype, v_ndims, v_dimids, &varid[varindex]); 
        ERR

        
        total_def_time += MPI_Wtime() - def_start_time;
        n_att = hdr_data->vars.value[i]->attrs.ndefined;
        assert(n_att == 0);
        // printf("\nn_att: %d\n", n_att);
        // for(k=0; k<n_att; k++){
        //     att_namelen = hdr_data->vars.value[i]->attrs.value[k]->name_len;
        //     att_xtype = hdr_data->vars.value[i]->attrs.value[k]->xtype;
        //     att_nelems = hdr_data->vars.value[i]->attrs.value[k]->nelems;
        //     err = ncmpi_put_att(ncid, varid[i], hdr_data->vars.value[i]->attrs.value[k]->name, att_xtype, 
        //     att_nelems, &hdr_data->vars.value[i]->attrs.value[k]->xvalue[0]); ERR
        // }
        free(v_dimids);

    }
    free(varid);
    return nerrs;
}




int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, nerrs=0;
    double end_time, start_time, start_time1, end_time1, start_time2, start_time3, max_time, min_time;
    double io_time, enddef_time, close_time, end_to_end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    struct hdr all_hdr;
    // struct hdr recv_hdr;
    // create_dummy_data(rank, &dummy);
    read_metadata_from_file(SOURCE_NAME, &all_hdr);
   // Print the created data for each process

    // printf("Total Header Size: %lld\n", all_hdr.xsz);
    // printf("Dimensions:\n");
    // for (int i = 0; i < 10; i++) {
    //     printf("  Name: %s, Size: %lld\n", all_hdr.dims.value[i]->name, all_hdr.dims.value[i]->size);
    // }

    // printf("Variables:\n");
    // for (int i = 0; i < 10; i++) {
    //     printf("Var %d  Name: %s, Type: %d, NumDims: %d", i, all_hdr.vars.value[i]->name,  all_hdr.vars.value[i]->xtype, 
    //     all_hdr.vars.value[i]->ndims);
    //     printf("    Dim IDs: ");
    //     for (int j = 0; j < all_hdr.vars.value[i]->ndims; j++) {
    //         printf("%d ", all_hdr.vars.value[i]->dimids[j]);
    //     }
    //     printf("\n");
    //     // printf("    Attributes:\n");
    //     // for (int k = 0; k < all_hdr.vars.value[i]->attrs.ndefined; k++) {
    //     //     printf("      Name: %s, Nelems: %lld, Type: %d\n", all_hdr.vars.value[i]->attrs.value[k]->name, 
    //     //     all_hdr.vars.value[i]->attrs.value[k]->nelems, all_hdr.vars.value[i]->attrs.value[k]->xtype);
    //     // }
    // }



    

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
    MPI_Info_set(info, "nc_hash_size_dim", "16777216");
    MPI_Info_set(info, "nc_hash_size_var", "8388608");

    snprintf(filename, sizeof(filename), "%.*s_%d.pnc", 
            (int)(sizeof(OUTPUT_NAME) - 5), OUTPUT_NAME, nproc);
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, info, &ncid); ERR
    double create_time = MPI_Wtime() - start_time;
    // MPI_Barrier(MPI_COMM_WORLD);

    // printf("rank %d, recv_displs: %d, recvcounts: %d \n",  rank, recv_displs[i], recvcounts[i]);
    int nvars = all_hdr.vars.ndefined;
    // int nvars = 10;
    int vars_per_process = nvars / nproc;
    int remainder = nvars % nproc;
    int start = rank * vars_per_process + (rank < remainder ? rank : remainder);
    int count = vars_per_process + (rank < remainder ? 1 : 0);
    // if (rank == 0) 
    //     printf("\ntotal var: %d", nvars);
    // printf("\nrank %d, start %d, count %d\n", rank, start, count);
    start_time1 = MPI_Wtime();
    define_hdr_nf(&all_hdr, ncid, rank, start, count);
    io_time = MPI_Wtime() - start_time1;
    
    
    // MPI_Barrier(MPI_COMM_WORLD);
    start_time2 = MPI_Wtime();
    err = ncmpi_enddef(ncid); ERR
    enddef_time = MPI_Wtime() - start_time2;

    // Clean up
    // MPI_Barrier(MPI_COMM_WORLD);
    start_time3 = MPI_Wtime();
    err = ncmpi_close(ncid); ERR
    end_time =  MPI_Wtime();
    close_time = end_time - start_time3;
    end_to_end_time = end_time - start_time;
    free_hdr(&all_hdr);


    double times[6] = {end_to_end_time, create_time,  io_time, enddef_time, close_time, total_def_time};
    char *names[6] = {"end-end", "create", "write", "enddef", "close", "def_dim/var"};
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