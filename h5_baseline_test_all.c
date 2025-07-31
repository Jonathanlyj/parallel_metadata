#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include <inttypes.h>
#include <sys/resource.h>
#include "baseline_ncx_app.h"




#define FAIL -1

#define DIRTY_BYTES_THRESHOLD (16 * 1024 * 1024) // 16 MB

const char *src_file = NULL;
const char *out_file = NULL;
double crt_start_time, total_crt_time=0;
int group_create_count = 0;
int dataset_create_count = 0;
//used to avoid open & closing same group for every dataset

char last_group_name[256] = ""; 
// Define a structure to represent a dataset
// typedef struct {
//     int name_len;
//     int dtype_size;
//     int dspace_size;
//     char *name;
//     char *dtype;
//     char *dspace;
// } h5_dataset;

// // Define a structure to represent a group
// typedef struct {
//     char *name;
//     int name_len;
//     int ndst;
//     h5_dataset **datasets;
// } h5_group;

// // Define a structure to represent an array of groups
// typedef struct {
//     int ngrps;
//     h5_group **groups;
// } h5_grouparray;


// void print_datatype(hid_t datatype_id) {
//     if (H5Tget_class(datatype_id) == H5T_INTEGER) {
//         printf("%"PRId64 "is Integer\n",datatype_id);
//     } else if (H5Tget_class(datatype_id) == H5T_FLOAT) {
//         printf("%"PRId64 "is Float\n",datatype_id);
//     } else if (H5Tget_class(datatype_id) == H5T_STRING) {
//         printf("%"PRId64 "is String\n",datatype_id);
//     } else {
//         printf("%"PRId64 "is Other\n",datatype_id);
//     }
// }
// Function to create a group

// Helper to extract group and dataset names from full variable name
void split_var_name(const char* full_name, char* group_name, char* dataset_name) {
    const char* ptr = full_name;
    int underscore_count = 0;
    const char* split_ptr = NULL;

    while (*ptr != '\0') {
        if (*ptr == '_') {
            underscore_count++;
            if (underscore_count == 5) {
                split_ptr = ptr;
                break;
            }
        }
        ptr++;
    }

    if (split_ptr) {
        int group_len = split_ptr - full_name;
        strncpy(group_name, full_name, group_len);
        group_name[group_len] = '\0';

        if (strncmp(group_name, "_grp", 4) == 0) {
            memmove(group_name, group_name + 4, strlen(group_name + 4) + 1);
        }

        strcpy(dataset_name, split_ptr + 1);
    } else {
        strcpy(group_name, "/");
        strcpy(dataset_name, full_name);
    }
}

static int deserialize_all_hdr(struct hdr **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (struct hdr *)malloc(sizeof(struct hdr));
        deserialize_hdr_meta(all_recv_hdr[i], all_collections_buffer + recv_displs[i], recvcounts[i]);
    }
    return 0;
}

// void new_group(const char *name, int name_len, int ndst, h5_group* group) {
//     group->name = strdup(name);
//     group->name_len = name_len;
//     group->ndst = ndst;
//     group->datasets = (h5_dataset**)malloc(ndst * sizeof(h5_dataset*));
//     return;
// }


// void new_dataset(const char *name, int name_len, char* dtype_tmp, size_t dtype_size, char* dspace_tmp, size_t dspace_size, h5_dataset* dataset){
//     dataset->name = strdup(name);
//     dataset->name_len = name_len;
//     dataset->dtype_size = dtype_size;
//     dataset->dspace_size = dspace_size;
//     dataset->dtype = malloc(dtype_size);
//     memcpy(dataset->dtype, dtype_tmp, dtype_size);
//     dataset->dspace = malloc(dspace_size);
//     memcpy(dataset->dspace, dspace_tmp, dspace_size);
// }

// // Function to free memory used by a group
// void free_dataset(h5_dataset *dataset) {
//     if (dataset == NULL) {
//         return;
//     }
//     free(dataset->name);
//     free(dataset->dtype);
//     free(dataset->dspace);
//     free(dataset);
// }
// void free_group(h5_group *group) {
//     free(group->name);
//     for (int i = 0; i < group->ndst; i++) {
//         free_dataset(group->datasets[i]);
//     }
//     free(group->datasets);
//     free(group);
// }

// // Function to free memory used by a group array
// void free_grouparray(h5_grouparray *local_meta) {
//     for (int i = 0; i < local_meta->ngrps; i++) {
//         free_group(local_meta->groups[i]);
//     }
//     free(local_meta->groups);
//     free(local_meta);
// }

/* ---------------------------------- Convert Metadata ----------------------------------------*/
// void read_metadata(h5_grouparray* local_meta){
//     hid_t       file_id, group_id, dst_id, datatype_id, dataspace_id, memspace;
//     herr_t      status;
//     hsize_t ngrp, ndst;
//     hsize_t *dims_out;
//     size_t dtype_size, dspace_size;
//     char *dtype_tmp, *dspace_tmp;
//    // Open the file
//     file_id = H5Fopen(src_file, H5F_ACC_RDONLY, H5P_DEFAULT);

//     // Get the number of groups in the root
    
//     H5Gget_num_objs(file_id, &ngrp);
    
//     // Initialize the group array
//     local_meta->ngrps = ngrp;
//     local_meta->groups = (h5_group**)malloc(ngrp * sizeof(h5_group*));

//     // Iterate over groups in the root
//     for (int i = 0; i < ngrp; i++) {
//         char group_name[256];
//         local_meta->groups[i] = (h5_group*)malloc(sizeof(h5_group));
        
//         H5Gget_objname_by_idx(file_id, (hsize_t)i, group_name, 256);

//         // Open the group
//         group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

//         // Create a group struct and add it to the group array

//         H5Gget_num_objs(group_id, &ndst);
//         new_group(group_name, strlen(group_name), ndst, local_meta->groups[i]);
//         for (int j = 0; j < ndst; j++) {
//             char dst_name[256];
//             local_meta->groups[i]->datasets[j] = malloc(sizeof(h5_dataset));
//             h5_dataset* dataset = local_meta->groups[i]->datasets[j];
//             H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);

//             dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
//             // Get dataspace and datatype_id
//             datatype_id = H5Dget_type(dst_id);
//             status = H5Tencode(datatype_id, NULL, &dtype_size);
//             dtype_tmp = malloc(dtype_size);
//             status = H5Tencode(datatype_id, dtype_tmp, &dtype_size);

//             dataspace_id = H5Dget_space(dst_id);
//             status = H5Sencode1(dataspace_id, NULL, &dspace_size);
//             dspace_tmp = malloc(dspace_size);
//             status = H5Sencode1(dataspace_id, dspace_tmp, &dspace_size);

//             new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
//             free(dtype_tmp);
//             free(dspace_tmp);

//             H5Tclose(datatype_id);
//             H5Sclose(dataspace_id);
//             H5Dclose(dst_id);
//         }
//         // Close the group
//         H5Gclose(group_id);
//     }
//     // Close/release resources
//     H5Fclose(file_id);
// }

// void read_metadata_parallel(h5_grouparray* local_meta, int rank, int nproc){
//     hid_t       file_id, plist_id, group_id, dst_id, datatype_id, dataspace_id, memspace;
//     herr_t      status;
//     hsize_t ngrp, ndst;
//     hsize_t *dims_out;
//     int grps_per_proc, remainder, start, count;
//     // plist_id = H5Pcreate(H5P_FILE_ACCESS);
//     // H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
//     size_t dtype_size, dspace_size;
//     char *dtype_tmp, *dspace_tmp;
//    // Open the file
//    // Note: need to use serial I/O version of because H5Gget_objname_by_idx in dset loop errors out
//     // file_id = H5Fopen(src_file, H5F_ACC_RDONLY, plist_id);
//     file_id = H5Fopen(src_file, H5F_ACC_RDONLY, H5P_DEFAULT);
//     assert(file_id != FAIL);
    
//     // Get the number of groups in the root
//     H5Gget_num_objs(file_id, &ngrp);
    
//     grps_per_proc = ngrp / nproc;
//     remainder = ngrp % nproc;
//     // Determine start and count based on rank
//     start = rank * grps_per_proc + (rank < remainder ? rank : remainder);
//     count = grps_per_proc + (rank < remainder ? 1 : 0);
//     // Initialize the group array
//     local_meta->ngrps = count;
//     local_meta->groups = (h5_group**)malloc(ngrp * sizeof(h5_group*));


//     // Iterate over groups in the root
//     for (int i = start; i < start + count; ++i) {
//         char group_name[256];

//         local_meta->groups[i - start] = (h5_group*)malloc(sizeof(h5_group));
//         H5Gget_objname_by_idx(file_id, (hsize_t)i, group_name, 256);
//         // Open the group
//         group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

//         // Create a group struct and add it to the group array

//         H5Gget_num_objs(group_id, &ndst);
//         new_group(group_name, strlen(group_name), ndst, local_meta->groups[i - start]);
//       for (int j = 0; j < ndst; j++) {
//             char dst_name[256];
//             h5_dataset* dataset;
//             // if (i >= start && i < start + count){
//             local_meta->groups[i - start]->datasets[j] = malloc(sizeof(h5_dataset));
//             dataset = local_meta->groups[i - start]->datasets[j];
//                 // }
//             H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);
//             dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
//             // Get dataspace and datatype_id
//             datatype_id = H5Dget_type(dst_id);
//             status = H5Tencode(datatype_id, NULL, &dtype_size);
//             dtype_tmp = malloc(dtype_size);
//             status = H5Tencode(datatype_id, dtype_tmp, &dtype_size);

//             dataspace_id = H5Dget_space(dst_id);
//             status = H5Sencode1(dataspace_id, NULL, &dspace_size);
//             dspace_tmp = malloc(dspace_size);
//             status = H5Sencode1(dataspace_id, dspace_tmp, &dspace_size);

//             // if (i >= start && i < start + count) new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
//             new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
//             free(dtype_tmp);
//             free(dspace_tmp);

//             H5Tclose(datatype_id);
//             H5Sclose(dataspace_id);
//             H5Dclose(dst_id);
//         }
//         // Close the group
//         H5Gclose(group_id);
//     }
//     // Close/release resources
//     // H5Pclose(plist_id);
//     H5Fclose(file_id);
// }



// /* ---------------------------------- Decode Metadata ----------------------------------------*/
// // Function to calculate the total size needed for the buffer
// size_t calculate_buffer_size(h5_grouparray *grouparray) {
//     size_t total_size = sizeof(int); // For ngrps
//     for (int i = 0; i < grouparray->ngrps; i++) {
//         h5_group *group = grouparray->groups[i];
//         total_size += sizeof(int); // For name_len
//         total_size += group->name_len; // For name
//         total_size += sizeof(int); // For ndst
//         for (int j = 0; j < group->ndst; j++) {
//             h5_dataset *dataset = group->datasets[j];
//             total_size += sizeof(int); // For name_len
//             total_size += dataset->name_len; // For name
//             total_size += sizeof(size_t); // For datatype
//             total_size += dataset->dtype_size; 
//             total_size += sizeof(size_t); // For dataspace
//             total_size += dataset->dspace_size; 

//         }
//     }
//     return total_size;
// }
// void serialize_dataset(h5_dataset *dataset, char **ptr) {
//     memcpy(*ptr, &(dataset->name_len), sizeof(int));
//     *ptr += sizeof(int);
//     memcpy(*ptr, dataset->name, dataset->name_len);
//     *ptr += dataset->name_len;
//     memcpy(*ptr, &(dataset->dtype_size), sizeof(size_t));
//     *ptr += sizeof(size_t);
//     memcpy(*ptr, dataset->dtype, dataset->dtype_size);
//     *ptr += dataset->dtype_size;
//     memcpy(*ptr, &(dataset->dspace_size), sizeof(size_t));
//     *ptr += sizeof(size_t);
//     memcpy(*ptr, dataset->dspace, dataset->dspace_size);
//     *ptr += dataset->dspace_size;
// }

// // Helper function to serialize a group
// void serialize_group(h5_group *group, char **ptr) {
//     memcpy(*ptr, &(group->name_len), sizeof(int));
//     *ptr += sizeof(int);
//     memcpy(*ptr, group->name, group->name_len);
//     *ptr += group->name_len;
//     memcpy(*ptr, &(group->ndst), sizeof(int));
//     *ptr += sizeof(int);

//     for (int j = 0; j < group->ndst; j++) {
//         serialize_dataset(group->datasets[j], ptr);
//     }
// }

// // Helper function to serialize a grouparray
// void serialize_grouparray(h5_grouparray *grouparray, void *buffer) {
//     char *ptr = (char *)buffer;
//     memcpy(ptr, &(grouparray->ngrps), sizeof(int));
//     ptr += sizeof(int);

//     for (int i = 0; i < grouparray->ngrps; i++) {
//         h5_group *group = grouparray->groups[i];
//         serialize_group(group, &ptr);
//     }
// }

// void deserialize_dataset(h5_dataset *dataset, char **ptr) {
//     memcpy(&(dataset->name_len), *ptr, sizeof(int));
//     *ptr += sizeof(int);
//     dataset->name = (char *)malloc(dataset->name_len + 1);
//     memcpy(dataset->name, *ptr, dataset->name_len);
//     dataset->name[dataset->name_len] = '\0';
//     *ptr += dataset->name_len;
//     memcpy(&(dataset->dtype_size), *ptr, sizeof(size_t));
//     *ptr += sizeof(size_t);
//     dataset->dtype = (char *)malloc(dataset->dtype_size);
//     memcpy(dataset->dtype, *ptr, dataset->dtype_size);
//     *ptr += dataset->dtype_size;
//     memcpy(&(dataset->dspace_size), *ptr, sizeof(size_t));
//     *ptr += sizeof(size_t);
//     dataset->dspace = (char *)malloc(dataset->dspace_size);
//     memcpy(dataset->dspace, *ptr, dataset->dspace_size);
//     *ptr += dataset->dspace_size;
// }

// // Helper function to deserialize a group
// void deserialize_group(h5_group *group, char **ptr) {
//     memcpy(&(group->name_len), *ptr, sizeof(int));
//     *ptr += sizeof(int);
//     group->name = (char *)malloc(group->name_len + 1);
//     memcpy(group->name, *ptr, group->name_len);
//     group->name[group->name_len] = '\0';
//     *ptr += group->name_len;
//     memcpy(&(group->ndst), *ptr, sizeof(int));
//     *ptr += sizeof(int);
//     group->datasets = (h5_dataset **)malloc(group->ndst * sizeof(h5_dataset *));
//     for (int j = 0; j < group->ndst; j++) {
//         group->datasets[j] = (h5_dataset *)malloc(sizeof(h5_dataset));
//         deserialize_dataset(group->datasets[j], ptr);
//     }
// }


// // Function to deserialize a grouparray from a buffer
// void deserialize_grouparray(h5_grouparray *grouparray, void *buffer){
//     char *ptr = (char *)buffer;
//     memcpy(&(grouparray->ngrps), ptr, sizeof(int));
//     ptr += sizeof(int);
//     grouparray->groups = (h5_group **)malloc(grouparray->ngrps * sizeof(h5_group *));
//     for (int i = 0; i < grouparray->ngrps; i++) {
//         grouparray->groups[i] = (h5_group *)malloc(sizeof(h5_group));
//         deserialize_group(grouparray->groups[i], &ptr);
//     }
// }



// static int deserialize_all_grouparray(h5_grouparray **all_recv_meta, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
//     for (int i=0; i< nproc; i++){
//         all_recv_meta[i]= (h5_grouparray *)malloc(sizeof(h5_grouparray));
//         deserialize_grouparray(all_recv_meta[i], all_collections_buffer + recv_displs[i]);
//         // Access group and dataset information
//     }
//     return 0;
// }

// static int free_all_grouparray(h5_grouparray **all_recv_meta, int nproc){
//     if (all_recv_meta != NULL){
//         for (int i=0; i< nproc; i++) free_grouparray(all_recv_meta[i]);
//         free(all_recv_meta);
//     }
//     return 0;
// }

int define_hdr_hdf5(struct hdr *hdr_data, hid_t file_id) {
    int nerrs = 0;
    int nvars = hdr_data->vars.ndefined;

    hid_t current_group_id = -1;
    // For tracking group context
    
    

    for (int i = 0; i < nvars; i++) {
        hdr_var *var = hdr_data->vars.value[i];
        char group_name[256], dataset_name[256];
        split_var_name(var->name, group_name, dataset_name);

        // Create or reuse group
        if (strcmp(group_name, last_group_name) != 0) {
            //need to create a new group, the group name has changed
            if (current_group_id >= 0)
                H5Gclose(current_group_id);

            current_group_id = H5Gcreate2(file_id, group_name,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            strcpy(last_group_name, group_name);
        } else if(i == 0){
            // If it's the first variable and the group is the same as the last one, we need to reopen the group
            current_group_id = H5Gopen2(file_id, group_name, H5P_DEFAULT);
        }
        // Map NC_ type to HDF5 native type
        hid_t h5type;
        switch (var->xtype) {
            case NC_INT:    h5type = H5T_NATIVE_INT; break;
            case NC_FLOAT:  h5type = H5T_NATIVE_FLOAT; break;
            case NC_DOUBLE: h5type = H5T_NATIVE_DOUBLE; break;
            case NC_CHAR:   h5type = H5T_NATIVE_CHAR; break;
            default:        h5type = H5T_NATIVE_INT; break;
        }

        // Create dataspace using dimension sizes
        hsize_t *dims = (hsize_t *)malloc(var->ndims * sizeof(hsize_t));
        for (int j = 0; j < var->ndims; j++) {
            dims[j] = hdr_data->dims.value[var->dimids[j]]->size;
        }
        hid_t space_id = H5Screate_simple(var->ndims, dims, NULL);
        free(dims);

        // Dataset creation property list
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_LATE);
        H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);


        


        // Create dataset in current group
        hid_t dset_id = H5Dcreate2(current_group_id, dataset_name, h5type, space_id,
                                   H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Pclose(dcpl_id);
        H5Sclose(space_id);

        

        // Create attributes
        for (int k = 0; k < var->attrs.ndefined; k++) {
            hdr_attr *att = var->attrs.value[k];

            hid_t att_type;
            switch (att->xtype) {
                case NC_INT:    att_type = H5T_NATIVE_INT; break;
                case NC_FLOAT:  att_type = H5T_NATIVE_FLOAT; break;
                case NC_DOUBLE: att_type = H5T_NATIVE_DOUBLE; break;
                case NC_CHAR:   att_type = H5T_C_S1; break;
                default:        att_type = H5T_NATIVE_INT; break;
            }

            hsize_t att_dims = att->nelems;
            hid_t att_space = H5Screate_simple(1, &att_dims, NULL);
            hid_t attr_id = H5Acreate2(dset_id, att->name, att_type, att_space,
                                       H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attr_id, att_type, att->xvalue);
            H5Sclose(att_space);
            H5Aclose(attr_id);
        }

        H5Dclose(dset_id);
    }
    H5Gclose(current_group_id);

    return nerrs;
}

static int free_all_hdr(struct hdr **all_recv_hdr, int nproc){
    if (all_recv_hdr != NULL){
        for (int i=0; i< nproc; i++) free_hdr_meta(all_recv_hdr[i]);
        free(all_recv_hdr);
    }
    return 0;
}




int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc, status, err, nerrs=0;
    hid_t outfile_id, plist_id, fcpl_id;
    double end_to_end_time, mpi_time, io_time, enddef_time, close_time, max_time, min_time;
    double start_time, start_time1, end_time1, end_time2, end_time3, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (argc < 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <source_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    src_file = argv[1];
    out_file = argv[2];

    if (rank == 0) {
        printf("Input file: %s\n", src_file);
        printf("Output file: %s\n", out_file);
    }
    struct hdr all_hdr;
    read_metadata_from_file(src_file, &all_hdr);
    struct hdr local_hdr;
    distribute_metadata(rank, nproc, &all_hdr, &local_hdr);
    free_hdr_meta(&all_hdr);
    // struct hdr recv_hdr;
    // create_local_hdr_data(rank, &local_hdr);

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
    status = serialize_hdr_meta(&local_hdr, send_buffer);


    


    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) malloc(nproc * sizeof(MPI_Offset));
    MPI_Allgather(&local_hdr.xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

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
    if(rank==0){
        printf("\nTotal buffer size: %f MB", total_recv_size_MB);
        printf("\nMax buffer size: %f MB", max_size_MB);
        printf("\nMin buffer size: %f MB \n", min_size_MB);
    }
    
    // printf("\nrank %d, local_hdr xsz %lld", rank, local_hdr.xsz);
    // Allocate buffer for receiving all header data
    char* all_collections_buffer = (char*) malloc(total_recv_size);
    int* recvcounts =  (int*)malloc(nproc * sizeof(int));
    for (int i = 0; i < nproc; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    // Phase 2: Communicate the actual header data
    // Before MPI_Allgatherv
    MPI_Allgatherv(send_buffer, local_hdr.xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);
    // Deserialize the received data and print if rank is 0
    
    int ncid, cmode;
    char filename[256];
    cmode = NC_64BIT_DATA | NC_CLOBBER;
    // size_t position = strlen(src_file) - 3;
    // strncpy(filename, src_file, position);
    // strcat(filename, "_new");
    // strcat(filename, src_file + position);
    // if (rank==0) printf("\n%s\n", out_file);
    struct hdr **all_recv_hdr = (struct hdr **)malloc(nproc * sizeof(struct hdr*));

    deserialize_all_hdr(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);
    
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    int block_size = 4 * 1024 * 1024;
    unsigned ik = 32;
    unsigned lk = 5;
    fcpl_id = H5Pcreate(H5P_FILE_CREATE);
    // H5Pset_istore_k(fcpl_id, 1024);
    H5Pset_sym_k(fcpl_id, ik, lk);
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    H5Pset_coll_metadata_write(plist_id, true);
    H5Pset_meta_block_size(plist_id, block_size);

    //Set dirty bytes threshold to 16 MB
    H5AC_cache_config_t mdc_config;
    // Initialize and retrieve current metadata cache configuration
    memset(&mdc_config, 0, sizeof(mdc_config));
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    H5Pget_mdc_config(plist_id, &mdc_config);

    // Set the dirty bytes threshold to 16 MB, so all processes will reach sync point until this amount of matadata cache is reached. Default: 256kb
    mdc_config.dirty_bytes_threshold = DIRTY_BYTES_THRESHOLD;
    // mdc_config.metadata_write_strategy = H5AC_METADATA_WRITE_STRATEGY__PROCESS_0_ONLY;

    // Apply the new configuration
    H5Pset_mdc_config(plist_id, &mdc_config);
    outfile_id = H5Fcreate(out_file, H5F_ACC_TRUNC, fcpl_id, plist_id);
    for (int i = 0; i < nproc; i++) {
        define_hdr_hdf5(all_recv_hdr[i], outfile_id);
    }
    
//     struct rusage r_usage;
//     getrusage(RUSAGE_SELF,&r_usage);
//   // Print the maximum resident set size used (in kilobytes).
//     printf("\n Memory usage: %ld kilobytes\n",r_usage.ru_maxrss);

    // create_metadata(new_meta, outfile_id);
    io_time = MPI_Wtime() - end_time1;
    H5Pclose(plist_id);
    H5Pclose(fcpl_id);
    free_all_hdr(all_recv_hdr, nproc);
    MPI_Barrier(MPI_COMM_WORLD);
    
    end_time3 = MPI_Wtime();
    H5Fclose(outfile_id);
    end_time = MPI_Wtime();
    close_time = end_time - end_time3;
    end_to_end_time = end_time - start_time;
    mpi_time = end_time1 - start_time1;
    double times[5] = {end_to_end_time, mpi_time, io_time, total_crt_time, close_time};
    char *names[5] = {"end-end", "mpi-phase", "write", "create", "close"};
    double max_times[5], min_times[5];


    MPI_Reduce(&times[0], &max_times[0], 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times[0], &min_times[0], 5, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) printf("ik: %u, lk: %u\n", ik, lk);
    for (int i = 0; i < 5; i++){
        if (rank == 0) {
            
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
    for (int i = 0; i < 5; i++){
        if (rank == 0) {
            printf("%f\n", names[i], max_times[i]);
            printf("%f\n", names[i], min_times[i]);
        }
    }
        if (rank == 0) {
        printf("H5Gcreate2 called %d times on rank 0\n", group_create_count);
        printf("H5Dcreate called %d times on rank 0\n", dataset_create_count);
    }
    // Free memory used by the group array
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    free(recvcounts);
    MPI_Finalize();
    return 0;
}
