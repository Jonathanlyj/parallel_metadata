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
#include "baseline_ncx.h" 



static int verbose;











    /* ---------------------------------- Metadata Generation ----------------------------------------*/

void create_dummy_data(int rank, struct hdr *dummy) {
    // Initialize xsz to the size of the hdr struct itself
    dummy->xsz = 0;
    int attrV_xsz, v_attrV_xsz;
    int status;

    // Create dummy dimensions
    dummy->dims.ndefined = rank;
    dummy->dims.value = (hdr_dim **)malloc(rank * sizeof(hdr_dim *));
    dummy->xsz += 2 * sizeof(uint32_t); // NC_Dimension and nelems
    for (int i = 0; i < rank; i++) {
        dummy->dims.value[i] = (hdr_dim *)malloc(sizeof(hdr_dim));
        dummy->dims.value[i]->size = i + 1;
        dummy->dims.value[i]->name_len = snprintf(NULL, 0, "dim_rank_%d_%d", rank, i);
        dummy->dims.value[i]->name = (char *)malloc(dummy->dims.value[i]->name_len + 1);
        sprintf(dummy->dims.value[i]->name, "dim_rank_%d_%d", rank, i);
        dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->dims.value[i]->name_len; // dim name
        dummy->xsz += sizeof(uint32_t); //size
    }

    // Create dummy global attributes
    dummy->attrs.ndefined = rank;
    dummy->attrs.value = (hdr_attr **)malloc(rank * sizeof(hdr_attr *));

    dummy->xsz += 2 * sizeof(uint32_t); // NC_Attribute and nelems
    
    for (int i = 0; i < rank; i++) {
        
        dummy->attrs.value[i] = (hdr_attr *)malloc(sizeof(hdr_attr));
        dummy->attrs.value[i]->nelems = i + 1;
        dummy->attrs.value[i]->xtype = NC_INT; // Using NC_INT for simplicity
        dummy->attrs.value[i]->name_len = snprintf(NULL, 0, "attr_rank_%d_%d", rank, i);
        dummy->attrs.value[i]->name = (char *)malloc(dummy->attrs.value[i]->name_len + 1);
        sprintf(dummy->attrs.value[i]->name, "attr_rank_%d_%d", rank, i);

        dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->attrs.value[i]->name_len; //attr name
        dummy->xsz += sizeof(uint32_t); // nc_type
        dummy->xsz += sizeof(uint32_t); // nelems
        status = xlen_nc_type(dummy->attrs.value[i]->xtype, &attrV_xsz);
        dummy->xsz += dummy->attrs.value[i]->nelems * attrV_xsz; // attr_value
    }

    // Create dummy variables
    dummy->vars.ndefined = rank; 
    dummy->vars.value = (hdr_var **)malloc(rank * sizeof(hdr_var *));
    dummy->xsz += 2 * sizeof(uint32_t); // NC_Variable and nelems
    for (int i = 0; i < rank; i++) {
        dummy->vars.value[i] = (hdr_var *)malloc(sizeof(hdr_var));
        dummy->vars.value[i]->xtype = NC_INT; // Using NC_INT for simplicity
        dummy->vars.value[i]->name_len = snprintf(NULL, 0, "var_rank_%d_%d", rank, i);
        dummy->vars.value[i]->name = (char *)malloc(dummy->vars.value[i]->name_len + 1);
        sprintf(dummy->vars.value[i]->name, "var_rank_%d_%d", rank, i);
        dummy->vars.value[i]->ndims = i + 1;
        dummy->vars.value[i]->dimids = (int *)malloc((i + 1) * sizeof(int));
        for (int j = 0; j <= i; j++) {
            dummy->vars.value[i]->dimids[j] = j % rank;
        }
        dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->vars.value[i]->name_len; //var name
        dummy->xsz += sizeof(uint32_t); //nelems of dim list
        dummy->xsz += sizeof(uint32_t) * dummy->vars.value[i]->ndims; // dimid list
        



        // Create dummy variable attributes
        dummy->vars.value[i]->attrs.ndefined = i;
        dummy->vars.value[i]->attrs.value = (hdr_attr **)malloc(i * sizeof(hdr_attr *));

        dummy->xsz += 2 * sizeof(uint32_t); // NC_Attribute and nelems

        for (int k = 0; k < i; k++) {
            dummy->vars.value[i]->attrs.value[k] = (hdr_attr *)malloc(sizeof(hdr_attr));
            dummy->vars.value[i]->attrs.value[k]->nelems = i + 1;
            dummy->vars.value[i]->attrs.value[k]->xtype = NC_INT; // Using NC_INT for simplicity
            dummy->vars.value[i]->attrs.value[k]->name_len = snprintf(NULL, 0, "var_attr_rank_%d_%d_%d", rank, i, k);
            dummy->vars.value[i]->attrs.value[k]->name = (char *)malloc(dummy->vars.value[i]->attrs.value[k]->name_len + 1);
            sprintf(dummy->vars.value[i]->attrs.value[k]->name, "var_attr_rank_%d_%d_%d", rank, i, k);

            dummy->xsz += sizeof(uint32_t) + sizeof(char) * dummy->vars.value[i]->attrs.value[k]->name_len; //attr name
            dummy->xsz += sizeof(uint32_t); // nc_type
            dummy->xsz += sizeof(uint32_t); // nelems
            status = xlen_nc_type(dummy->vars.value[i]->attrs.value[k]->xtype, &v_attrV_xsz);
            dummy->xsz += dummy->vars.value[i]->attrs.value[k]->nelems * attrV_xsz; // attr_value

        }
    }
}



int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    struct hdr dummy;
    create_dummy_data(rank, &dummy);

   // Print the created data for each process
    printf("Rank %d:\n", rank);
    printf("Total Header Size: %lld\n", dummy.xsz);
    printf("Dimensions:\n");
    for (int i = 0; i < dummy.dims.ndefined; i++) {
        printf("  Name: %s, Size: %lld\n", dummy.dims.value[i]->name, dummy.dims.value[i]->size);
    }

    printf("Global Attributes:\n");
    for (int i = 0; i < dummy.attrs.ndefined; i++) {
        printf("  Name: %s, Nelems: %lld, Type: %d\n", dummy.attrs.value[i]->name, dummy.attrs.value[i]->nelems, dummy.attrs.value[i]->xtype);
    }

    printf("Variables:\n");
    for (int i = 0; i < dummy.vars.ndefined; i++) {
        printf("  Name: %s, Type: %d, NumDims: %d\n", dummy.vars.value[i]->name,  dummy.vars.value[i]->xtype, dummy.vars.value[i]->ndims);
        printf("    Dim IDs: ");
        for (int j = 0; j < dummy.vars.value[i]->ndims; j++) {
            printf("%d ", dummy.vars.value[i]->dimids[j]);
        }
        printf("\n");

        printf("    Attributes:\n");
        for (int k = 0; k < dummy.vars.value[i]->attrs.ndefined; k++) {
            printf("      Name: %s, Nelems: %lld, Type: %d\n", dummy.vars.value[i]->attrs.value[k]->name, dummy.vars.value[i]->attrs.value[k]->nelems, dummy.vars.value[i]->attrs.value[k]->xtype);
        }
    }

    MPI_Finalize();

    // // Example: Each process has a number of dimensions equal to its rank + 1
    // int num_dimensions = rank + 1;

    // // Allocate and fill dimension collection for the current process
    // DimensionCollection send_collection;
    // send_collection.num_dimensions = num_dimensions;
    // send_collection.dimension_sizes = (int*) malloc(num_dimensions * sizeof(int));
    // send_collection.name_lengths = (int*) malloc(num_dimensions * sizeof(int));
    // send_collection.total_names_length = 0;
    // char* names_concat = (char*) malloc(num_dimensions * 50 * sizeof(char)); // Each name up to 50 chars
    // names_concat[0] = '\0';
    // for (int i = 0; i < num_dimensions; ++i) {
    //     char dimension_name[50];
    //     sprintf(dimension_name, "Process_%d_Dim_%d", rank, i);
    //     send_collection.dimension_sizes[i] = rank * 10 + i;
    //     send_collection.name_lengths[i] = strlen(dimension_name) + 1;
    //     send_collection.total_names_length += send_collection.name_lengths[i];
    //     strcat(names_concat, dimension_name);
    //     strcat(names_concat, ";"); // Separator
    // }
    // send_collection.dimension_names = names_concat;

    // // Serialize the send_collection
    // int my_collection_size = sizeof(int) * (2 + 2 * send_collection.num_dimensions) + sizeof(char) * send_collection.total_names_length;
    // char* send_buffer = (char*) malloc(my_collection_size);
    // serialize(&send_collection, send_buffer);

    // // Phase 1: Communicate the sizes of the DimensionCollection structure for each process
    // int* all_collection_sizes = (int*) malloc(nprocs * sizeof(int));
    // MPI_Allgather(&my_collection_size, 1, MPI_INT, all_collection_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    // // Calculate displacements for the second phase
    // int* recv_displs = (int*) malloc(nprocs * sizeof(int));
    // int total_recv_size = all_collection_sizes[0];
    // recv_displs[0] = 0;
    // for (int i = 1; i < nprocs; ++i) {
    //     recv_displs[i] = recv_displs[i - 1] + all_collection_sizes[i - 1];
    //     total_recv_size += all_collection_sizes[i];
    // }

    // // Allocate buffer for receiving all DimensionCollection data
    // char* all_collections_buffer = (char*) malloc(total_recv_size);

    // // Phase 2: Communicate the actual DimensionCollection data
    // MPI_Allgatherv(send_buffer, my_collection_size, MPI_BYTE,
    //                all_collections_buffer, all_collection_sizes, recv_displs, MPI_BYTE, MPI_COMM_WORLD);

    // printf("The collected buffer is: %s\n", all_collections_buffer);
    // Deserialize the received data
    // For example, deserialize the data for the first process (process 0)

// Deserialize the received data and print if rank is 0
// if (rank == 0){
//     for (int i = 0; i < nprocs; ++i) {
//         DimensionCollection received_collection;
//         deserialize(all_collections_buffer + recv_displs[i], &received_collection);

//         printf("Process %d:\n", i);
//         printf("  Total Name Length: %d\n", received_collection.total_names_length);
//         printf("  Number of Dimensions: %d\n", received_collection.num_dimensions);
//         char* name_ptr = received_collection.dimension_names;
//         for (int j = 0; j < received_collection.num_dimensions; ++j) {
//             printf("    Dimension %s: Size = %d, Name Length = %d\n",
//                    name_ptr, received_collection.dimension_sizes[j], received_collection.name_lengths[j]);
//             name_ptr += received_collection.name_lengths[j]; 
//         }
        
//         // Free the arrays inside received_collection
//         free(received_collection.dimension_sizes);
//         free(received_collection.name_lengths);
//         free(received_collection.dimension_names);
//     }
// }
    // ... Process the deserialized data ...

    // Clean up
    // free(send_buffer);
    // free(all_collections_buffer);
    // free(all_collection_sizes);
    // free(recv_displs);
    // free(send_collection.dimension_sizes);
    // free(send_collection.name_lengths);
    // free(send_collection.dimension_names);
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