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
#include "baseline_ncx_app.h" 



static int verbose;


#define ERR {if(err!=NC_NOERR){printf("Error at %s:%d : %s\n", __FILE__,__LINE__, ncmpi_strerror(err));nerrs++;}}


    /* ---------------------------------- Metadata Generation ----------------------------------------*/


int define_hdr(struct hdr *hdr_data, int ncid){
    //define dimensions
    int ndims= hdr_data->dims.ndefined;
    int *dimid = (int *)malloc(ndims * sizeof(int));
    int i,j,k,nerrs=0;
    int err;

     for (i=0; i<ndims; i++){
        err = ncmpi_def_dim(ncid, hdr_data->dims.value[i]->name,  hdr_data->dims.value[i]->size, &dimid[i]); ERR
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
        err = ncmpi_def_var(ncid, hdr_data->vars.value[i]->name, xtype, v_ndims,  v_dimids, &varid[i]); ERR
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


void create_dummy_data(int rank, struct hdr *dummy) {
    // Initialize xsz to the size of the hdr struct itself
    dummy->xsz = 0;
    int attrV_xsz, v_attrV_xsz;
    int status;
    int att_nelems;
    int* att_array;

    // Create dummy dimension
    dummy->dims.ndefined = 2;
    dummy->dims.value = (hdr_dim **)malloc(rank * sizeof(hdr_dim *));
    dummy->xsz += 2 * sizeof(uint32_t); // NC_Dimension and nelems
    for (int i = 0; i < dummy->dims.ndefined; i++) {
        dummy->dims.value[i] = (hdr_dim *)malloc(sizeof(hdr_dim));
        dummy->dims.value[i]->size = i + 1;
        dummy->dims.value[i]->name_len = snprintf(NULL, 0, "dim_rank_%d_%d", rank, i);
        dummy->dims.value[i]->name = (char *)malloc(dummy->dims.value[i]->name_len + 1);
        sprintf(dummy->dims.value[i]->name, "dim_rank_%d_%d", rank, i);
        dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->dims.value[i]->name_len; // dim name
        dummy->xsz += sizeof(uint32_t); //size
        
    }


    // Create dummy variables
    dummy->vars.ndefined = 2; 
    dummy->vars.value = (hdr_var **)malloc(dummy->vars.ndefined * sizeof(hdr_var *));
    dummy->xsz += 2 * sizeof(uint32_t); // NC_Variable and ndefined
    for (int i = 0; i < dummy->vars.ndefined; i++) {
        dummy->vars.value[i] = (hdr_var *)malloc(sizeof(hdr_var));
        dummy->vars.value[i]->xtype = NC_INT; // Using NC_INT for simplicity
        dummy->vars.value[i]->name_len = snprintf(NULL, 0, "var_rank_%d_%d", rank, i);
        dummy->vars.value[i]->name = (char *)malloc(dummy->vars.value[i]->name_len + 1);
        sprintf(dummy->vars.value[i]->name, "var_rank_%d_%d", rank, i);
        dummy->vars.value[i]->ndims = 2;
        dummy->vars.value[i]->dimids = (int *)malloc((i + 1) * sizeof(int));
        for (int j = 0; j <= dummy->vars.value[i]->ndims; j++) {
            dummy->vars.value[i]->dimids[j] = j;
        }
        dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->vars.value[i]->name_len; //var name
        dummy->xsz += sizeof(uint32_t); //xtype
        dummy->xsz += sizeof(uint32_t); //nelems of dim list
        dummy->xsz += sizeof(uint32_t) * dummy->vars.value[i]->ndims; // dimid list

        



    //     // Create dummy variable attributes
        dummy->vars.value[i]->attrs.ndefined = i;
        dummy->vars.value[i]->attrs.value = (hdr_attr **)malloc(i * sizeof(hdr_attr *));
        dummy->xsz += 2 * sizeof(uint32_t); // NC_Attribute and ndefined

        for (int k = 0; k < i; k++) {
            dummy->vars.value[i]->attrs.value[k] = (hdr_attr *)malloc(sizeof(hdr_attr));
            dummy->vars.value[i]->attrs.value[k]->nelems = i + 1;
            dummy->vars.value[i]->attrs.value[k]->xtype = NC_INT; // Using NC_INT for simplicity
            dummy->vars.value[i]->attrs.value[k]->name_len = snprintf(NULL, 0, "var_attr_rank_%d_%d_%d", rank, i, k);
            dummy->vars.value[i]->attrs.value[k]->name = (char *)malloc(dummy->vars.value[i]->attrs.value[k]->name_len + 1);
            sprintf(dummy->vars.value[i]->attrs.value[k]->name, "var_attr_rank_%d_%d_%d", rank, i, k);
            att_nelems = dummy->vars.value[i]->attrs.value[k]->nelems;
            att_array = (int*)malloc(att_nelems * sizeof(int));
            for (int m = 0; m < att_nelems; m++) {
                    att_array[m] = rank;
                }
            dummy->vars.value[i]->attrs.value[k]->xvalue = att_array;


            dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->vars.value[i]->attrs.value[k]->name_len; //attr name
            dummy->xsz += sizeof(uint32_t); // nc_type
            dummy->xsz += sizeof(uint32_t); // nelems
            status = xlen_nc_type(dummy->vars.value[i]->attrs.value[k]->xtype, &v_attrV_xsz);
            dummy->xsz += dummy->vars.value[i]->attrs.value[k]->nelems * v_attrV_xsz; // attr_value

        }
    }


}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size, status, err, nerrs=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    struct hdr dummy;
    struct hdr recv_hdr;
    create_dummy_data(rank, &dummy);
//    // Print the created data for each process
//     printf("\nRank %d:\n", rank);
//     printf("Total Header Size: %lld\n", dummy.xsz);
//     printf("Dimensions:\n");
//     for (int i = 0; i < dummy.dims.ndefined; i++) {
//         printf("  Name: %s, Size: %lld\n", dummy.dims.value[i]->name, dummy.dims.value[i]->size);
//     }

//     printf("Variables:\n");
//     for (int i = 0; i < dummy.vars.ndefined; i++) {
//         printf("  Name: %s, Type: %d, NumDims: %d\n", dummy.vars.value[i]->name,  dummy.vars.value[i]->xtype, 
//         dummy.vars.value[i]->ndims);
//         printf("    Dim IDs: ");
//         for (int j = 0; j < dummy.vars.value[i]->ndims; j++) {
//             printf("%d ", dummy.vars.value[i]->dimids[j]);
//         }
//         printf("\n");
//         printf("    Attributes:\n");
//         for (int k = 0; k < dummy.vars.value[i]->attrs.ndefined; k++) {
//             printf("      Name: %s, Nelems: %lld, Type: %d\n", dummy.vars.value[i]->attrs.value[k]->name, 
//             dummy.vars.value[i]->attrs.value[k]->nelems, dummy.vars.value[i]->attrs.value[k]->xtype);
//         }
//     }
//     printf("rank %d, buffer size: %lld \n", rank, dummy.xsz);
    char* send_buffer = (char*) malloc(dummy.xsz);
    status = serialize_hdr(&dummy, send_buffer);
    // printf("rank %d, buffer: %s", rank, send_buffer);




    


    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) malloc(size * sizeof(MPI_Offset));
    MPI_Allgather(&dummy.xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

    // Calculate displacements for the second phase
    int* recv_displs = (int*) malloc(size * sizeof(int));
    int total_recv_size = all_collection_sizes[0];
    recv_displs[0] = 0;

    
    for (int i = 1; i < size; ++i) {
        recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
        total_recv_size += all_collection_sizes[i];
        
    }
    
    // printf("\nrank %d, dummy xsz %lld", rank, dummy.xsz);
    // Allocate buffer for receiving all header data
    char* all_collections_buffer = (char*) malloc(total_recv_size);

    int* recvcounts =  (int*)malloc(size * sizeof(int));

    for (int i = 0; i < size; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    // Phase 2: Communicate the actual header data
    // Before MPI_Allgatherv
    MPI_Allgatherv(send_buffer, dummy.xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);



    // Deserialize the received data
    // For example, deserialize the data for the first process (process 0)

// Deserialize the received data and print if rank is 0
    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_OFFSET| NC_CLOBBER;
    strcpy(filename, "testfile.nc");
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); ERR
    for (int i = 0; i < size; ++i) {
        struct hdr recv_hdr;
        // printf("rank %d, recv_displs: %d, recvcounts: %d \n",  rank, recv_displs[i], recvcounts[i]);
        deserialize_hdr(&recv_hdr, all_collections_buffer + recv_displs[i], recvcounts[i]);
        define_hdr(&recv_hdr, ncid);
    }


    // Clean up
    free(send_buffer);
    // free(all_collections_buffer);
    free(all_collection_sizes);
    // free(recv_displs);
    // free(send_collection.dimension_sizes);
    // free(send_collection.name_lengths);
    // free(send_collection.dimension_names);
    // ... Also, free the arrays inside received_collection ...
    err = ncmpi_close(ncid); ERR
    MPI_Finalize();
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
    // for (i=0; i<size; i++){
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