#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include <inttypes.h>


#define SRC_FILE "nue_slice_panoptic_hdf_merged_meta.h5"
// #define SRC_FILE "/files2/scratch/yll6162/parallel_metadata/script/test/nue_slice_graphs.0001.h5"
// #define SRC_FILE "h5_example.h5"
#define OUT_FILE "h5_baseline_test_all.h5"
#define FAIL -1
double crt_start_time, total_crt_time=0;
// Define a structure to represent a dataset
typedef struct {
    int name_len;
    int dtype_size;
    int dspace_size;
    char *name;
    char *dtype;
    char *dspace;
} h5_dataset;

// Define a structure to represent a group
typedef struct {
    char *name;
    int name_len;
    int ndst;
    h5_dataset **datasets;
} h5_group;

// Define a structure to represent an array of groups
typedef struct {
    int ngrps;
    h5_group **groups;
} h5_grouparray;


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
void new_group(const char *name, int name_len, int ndst, h5_group* group) {
    group->name = strdup(name);
    group->name_len = name_len;
    group->ndst = ndst;
    group->datasets = (h5_dataset**)malloc(ndst * sizeof(h5_dataset*));
    return;
}


void new_dataset(const char *name, int name_len, char* dtype_tmp, size_t dtype_size, char* dspace_tmp, size_t dspace_size, h5_dataset* dataset){
    dataset->name = strdup(name);
    dataset->name_len = name_len;
    dataset->dtype_size = dtype_size;
    dataset->dspace_size = dspace_size;
    dataset->dtype = malloc(dtype_size);
    memcpy(dataset->dtype, dtype_tmp, dtype_size);
    dataset->dspace = malloc(dspace_size);
    memcpy(dataset->dspace, dspace_tmp, dspace_size);
}

// Function to free memory used by a group
void free_dataset(h5_dataset *dataset) {
    if (dataset == NULL) {
        return;
    }
    free(dataset->name);
    free(dataset->dtype);
    free(dataset->dspace);
    free(dataset);
}
void free_group(h5_group *group) {
    free(group->name);
    for (int i = 0; i < group->ndst; i++) {
        free_dataset(group->datasets[i]);
    }
    free(group->datasets);
    free(group);
}

// Function to free memory used by a group array
void free_grouparray(h5_grouparray *local_meta) {
    for (int i = 0; i < local_meta->ngrps; i++) {
        free_group(local_meta->groups[i]);
    }
    free(local_meta->groups);
    free(local_meta);
}

/* ---------------------------------- Read Metadata ----------------------------------------*/
void read_metadata(h5_grouparray* local_meta){
    hid_t       file_id, group_id, dst_id, datatype_id, dataspace_id, memspace;
    herr_t      status;
    hsize_t ngrp, ndst;
    hsize_t *dims_out;
    size_t dtype_size, dspace_size;
    char *dtype_tmp, *dspace_tmp;
   // Open the file
    file_id = H5Fopen(SRC_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Get the number of groups in the root
    
    H5Gget_num_objs(file_id, &ngrp);
    
    // Initialize the group array
    local_meta->ngrps = ngrp;
    local_meta->groups = (h5_group**)malloc(ngrp * sizeof(h5_group*));

    // Iterate over groups in the root
    for (int i = 0; i < ngrp; i++) {
        char group_name[256];
        local_meta->groups[i] = (h5_group*)malloc(sizeof(h5_group));
        
        H5Gget_objname_by_idx(file_id, (hsize_t)i, group_name, 256);

        // Open the group
        group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

        // Create a group struct and add it to the group array

        H5Gget_num_objs(group_id, &ndst);
        new_group(group_name, strlen(group_name), ndst, local_meta->groups[i]);
        for (int j = 0; j < ndst; j++) {
            char dst_name[256];
            local_meta->groups[i]->datasets[j] = malloc(sizeof(h5_dataset));
            h5_dataset* dataset = local_meta->groups[i]->datasets[j];
            H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);

            dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
            // Get dataspace and datatype_id
            datatype_id = H5Dget_type(dst_id);
            status = H5Tencode(datatype_id, NULL, &dtype_size);
            dtype_tmp = malloc(dtype_size);
            status = H5Tencode(datatype_id, dtype_tmp, &dtype_size);

            dataspace_id = H5Dget_space(dst_id);
            status = H5Sencode1(dataspace_id, NULL, &dspace_size);
            dspace_tmp = malloc(dspace_size);
            status = H5Sencode1(dataspace_id, dspace_tmp, &dspace_size);

            new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
            free(dtype_tmp);
            free(dspace_tmp);

            H5Tclose(datatype_id);
            H5Sclose(dataspace_id);
            H5Dclose(dst_id);
        }
        // Close the group
        H5Gclose(group_id);
    }
    // Close/release resources
    H5Fclose(file_id);
}

void read_metadata_parallel(h5_grouparray* local_meta, int rank, int nproc){
    hid_t       file_id, plist_id, group_id, dst_id, datatype_id, dataspace_id, memspace;
    herr_t      status;
    hsize_t ngrp, ndst;
    hsize_t *dims_out;
    int grps_per_proc, remainder, start, count;
    // plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    size_t dtype_size, dspace_size;
    char *dtype_tmp, *dspace_tmp;
   // Open the file
   // Note: need to use serial I/O version of because H5Gget_objname_by_idx in dset loop errors out
    // file_id = H5Fopen(SRC_FILE, H5F_ACC_RDONLY, plist_id);
    file_id = H5Fopen(SRC_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file_id != FAIL);
    
    // Get the number of groups in the root
    H5Gget_num_objs(file_id, &ngrp);
    
    grps_per_proc = ngrp / nproc;
    remainder = ngrp % nproc;
    // Determine start and count based on rank
    start = rank * grps_per_proc + (rank < remainder ? rank : remainder);
    count = grps_per_proc + (rank < remainder ? 1 : 0);
    // Initialize the group array
    local_meta->ngrps = count;
    local_meta->groups = (h5_group**)malloc(ngrp * sizeof(h5_group*));


    // Iterate over groups in the root
    for (int i = start; i < start + count; ++i) {
        char group_name[256];

        local_meta->groups[i - start] = (h5_group*)malloc(sizeof(h5_group));
        H5Gget_objname_by_idx(file_id, (hsize_t)i, group_name, 256);
        // Open the group
        group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);

        // Create a group struct and add it to the group array

        H5Gget_num_objs(group_id, &ndst);
        new_group(group_name, strlen(group_name), ndst, local_meta->groups[i - start]);
      for (int j = 0; j < ndst; j++) {
            char dst_name[256];
            h5_dataset* dataset;
            // if (i >= start && i < start + count){
            local_meta->groups[i - start]->datasets[j] = malloc(sizeof(h5_dataset));
            dataset = local_meta->groups[i - start]->datasets[j];
                // }
            H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);
            dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
            // Get dataspace and datatype_id
            datatype_id = H5Dget_type(dst_id);
            status = H5Tencode(datatype_id, NULL, &dtype_size);
            dtype_tmp = malloc(dtype_size);
            status = H5Tencode(datatype_id, dtype_tmp, &dtype_size);

            dataspace_id = H5Dget_space(dst_id);
            status = H5Sencode1(dataspace_id, NULL, &dspace_size);
            dspace_tmp = malloc(dspace_size);
            status = H5Sencode1(dataspace_id, dspace_tmp, &dspace_size);

            // if (i >= start && i < start + count) new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
            new_dataset(dst_name, strlen(dst_name),  dtype_tmp, dtype_size, dspace_tmp, dspace_size, dataset);
            free(dtype_tmp);
            free(dspace_tmp);

            H5Tclose(datatype_id);
            H5Sclose(dataspace_id);
            H5Dclose(dst_id);
        }
        // Close the group
        H5Gclose(group_id);
    }
    // Close/release resources
    // H5Pclose(plist_id);
    H5Fclose(file_id);
}

/* ---------------------------------- Decode Metadata ----------------------------------------*/
// Function to calculate the total size needed for the buffer
size_t calculate_buffer_size(h5_grouparray *grouparray) {
    size_t total_size = sizeof(int); // For ngrps
    for (int i = 0; i < grouparray->ngrps; i++) {
        h5_group *group = grouparray->groups[i];
        total_size += sizeof(int); // For name_len
        total_size += group->name_len; // For name
        total_size += sizeof(int); // For ndst
        for (int j = 0; j < group->ndst; j++) {
            h5_dataset *dataset = group->datasets[j];
            total_size += sizeof(int); // For name_len
            total_size += dataset->name_len; // For name
            total_size += sizeof(size_t); // For datatype
            total_size += dataset->dtype_size; 
            total_size += sizeof(size_t); // For dataspace
            total_size += dataset->dspace_size; 

        }
    }
    return total_size;
}
void serialize_dataset(h5_dataset *dataset, char **ptr) {
    memcpy(*ptr, &(dataset->name_len), sizeof(int));
    *ptr += sizeof(int);
    memcpy(*ptr, dataset->name, dataset->name_len);
    *ptr += dataset->name_len;
    memcpy(*ptr, &(dataset->dtype_size), sizeof(size_t));
    *ptr += sizeof(size_t);
    memcpy(*ptr, dataset->dtype, dataset->dtype_size);
    *ptr += dataset->dtype_size;
    memcpy(*ptr, &(dataset->dspace_size), sizeof(size_t));
    *ptr += sizeof(size_t);
    memcpy(*ptr, dataset->dspace, dataset->dspace_size);
    *ptr += dataset->dspace_size;
}

// Helper function to serialize a group
void serialize_group(h5_group *group, char **ptr) {
    memcpy(*ptr, &(group->name_len), sizeof(int));
    *ptr += sizeof(int);
    memcpy(*ptr, group->name, group->name_len);
    *ptr += group->name_len;
    memcpy(*ptr, &(group->ndst), sizeof(int));
    *ptr += sizeof(int);

    for (int j = 0; j < group->ndst; j++) {
        serialize_dataset(group->datasets[j], ptr);
    }
}

// Helper function to serialize a grouparray
void serialize_grouparray(h5_grouparray *grouparray, void *buffer) {
    char *ptr = (char *)buffer;
    memcpy(ptr, &(grouparray->ngrps), sizeof(int));
    ptr += sizeof(int);

    for (int i = 0; i < grouparray->ngrps; i++) {
        h5_group *group = grouparray->groups[i];
        serialize_group(group, &ptr);
    }
}

void deserialize_dataset(h5_dataset *dataset, char **ptr) {
    memcpy(&(dataset->name_len), *ptr, sizeof(int));
    *ptr += sizeof(int);
    dataset->name = (char *)malloc(dataset->name_len + 1);
    memcpy(dataset->name, *ptr, dataset->name_len);
    dataset->name[dataset->name_len] = '\0';
    *ptr += dataset->name_len;
    memcpy(&(dataset->dtype_size), *ptr, sizeof(size_t));
    *ptr += sizeof(size_t);
    dataset->dtype = (char *)malloc(dataset->dtype_size);
    memcpy(dataset->dtype, *ptr, dataset->dtype_size);
    *ptr += dataset->dtype_size;
    memcpy(&(dataset->dspace_size), *ptr, sizeof(size_t));
    *ptr += sizeof(size_t);
    dataset->dspace = (char *)malloc(dataset->dspace_size);
    memcpy(dataset->dspace, *ptr, dataset->dspace_size);
    *ptr += dataset->dspace_size;
}

// Helper function to deserialize a group
void deserialize_group(h5_group *group, char **ptr) {
    memcpy(&(group->name_len), *ptr, sizeof(int));
    *ptr += sizeof(int);
    group->name = (char *)malloc(group->name_len + 1);
    memcpy(group->name, *ptr, group->name_len);
    group->name[group->name_len] = '\0';
    *ptr += group->name_len;
    memcpy(&(group->ndst), *ptr, sizeof(int));
    *ptr += sizeof(int);
    group->datasets = (h5_dataset **)malloc(group->ndst * sizeof(h5_dataset *));
    for (int j = 0; j < group->ndst; j++) {
        group->datasets[j] = (h5_dataset *)malloc(sizeof(h5_dataset));
        deserialize_dataset(group->datasets[j], ptr);
    }
}


// Function to deserialize a grouparray from a buffer
void deserialize_grouparray(h5_grouparray *grouparray, void *buffer){
    char *ptr = (char *)buffer;
    memcpy(&(grouparray->ngrps), ptr, sizeof(int));
    ptr += sizeof(int);
    grouparray->groups = (h5_group **)malloc(grouparray->ngrps * sizeof(h5_group *));
    for (int i = 0; i < grouparray->ngrps; i++) {
        grouparray->groups[i] = (h5_group *)malloc(sizeof(h5_group));
        deserialize_group(grouparray->groups[i], &ptr);
    }
}



static int deserialize_all_grouparray(h5_grouparray **all_recv_meta, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_meta[i]= (h5_grouparray *)malloc(sizeof(h5_grouparray));
        deserialize_grouparray(all_recv_meta[i], all_collections_buffer + recv_displs[i]);
        // Access group and dataset information
    }
    return 0;
}

static int free_all_grouparray(h5_grouparray **all_recv_meta, int nproc){
    if (all_recv_meta != NULL){
        for (int i=0; i< nproc; i++) free_grouparray(all_recv_meta[i]);
        free(all_recv_meta);
    }
    return 0;
}

/* ---------------------------------- Write Metadata ----------------------------------------*/
void create_metadata(h5_grouparray* local_meta, hid_t file_id) {
    hid_t group_id, dataset_id, dataspace_id, datatype_id, dcpl_id;
    H5D_space_status_t space_status;
    hsize_t            storage_size;
    herr_t status;
    // Iterate over each group in the local_meta
    for (int i = 0; i < local_meta->ngrps; i++) {
        h5_group* group = local_meta->groups[i];
        // Create a group in the HDF5 file
        // crt_start_time = MPI_Wtime();
        // group_id = H5Gcreate2(file_id, group->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // total_crt_time += MPI_Wtime() - crt_start_time;

        // Iterate over each dataset in the group
        for (int j = 0; j < group->ndst; j++) {
            h5_dataset* dataset = group->datasets[j];
            char* flat_dataset_name = (char*)malloc(strlen(group->name) + strlen(dataset->name) + 2);
            strcpy(flat_dataset_name, group->name);
            strcat(flat_dataset_name, dataset->name);
            // Create dataspace
            dataspace_id = H5Sdecode(dataset->dspace);
            // Create datatypes
            datatype_id = H5Tdecode(dataset->dtype);
            // Create dataset
            dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            // Set allocation time to late
            status = H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_LATE);//no effect under parallel-io mode
            // // Turn off auto-fill for the dataset
            status = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
            crt_start_time = MPI_Wtime();
            dataset_id = H5Dcreate(file_id, flat_dataset_name, datatype_id, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            total_crt_time += MPI_Wtime() - crt_start_time;

            // check storage space allocated
            // status       = H5Dget_space_status(dataset_id, &space_status);
            // storage_size = H5Dget_storage_size(dataset_id);
            // printf("Space for %s has%sbeen allocated.\n", dataset->name,
            //     space_status == H5D_SPACE_STATUS_ALLOCATED ? " " : " not ");
            // printf("Storage size for %s is: %ld bytes.\n", dataset->name, (long)storage_size);
                        // Close resources
            H5Dclose(dataset_id);
            H5Pclose(dcpl_id);
            H5Tclose(datatype_id);
            H5Sclose(dataspace_id);
        }

        // Close the group
        // H5Gclose(group_id);
    }
}

int create_all_metadata(h5_grouparray **all_recv_meta, int nproc, hid_t file_id){
    for (int i=0; i<nproc; i++){
        h5_grouparray *hdr_data = all_recv_meta[i];
        create_metadata(hdr_data, file_id);
    }
    return 0;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, nproc;
    hid_t outfile_id, plist_id ;
    h5_grouparray* local_meta;
    size_t local_xsz, global_xsz;
    hsize_t block_size;
    
    local_meta = (h5_grouparray*)malloc(sizeof(h5_grouparray));
    
    char *send_buffer;

    double end_to_end_time, mpi_time, io_time, enddef_time, close_time, max_time, min_time;
    double start_time, start_time1, end_time1, end_time2, end_time3, end_time;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
 
    read_metadata_parallel(local_meta, rank, nproc);
    // read_metadata(local_meta);

    MPI_Barrier(MPI_COMM_WORLD);
    start_time = start_time1 = MPI_Wtime();
    local_xsz = calculate_buffer_size(local_meta);

    send_buffer = (char*) malloc(local_xsz);

    serialize_grouparray(local_meta, send_buffer);

    
    // Phase 1: Communicate the sizes of the header structure for each process
    MPI_Offset* all_collection_sizes = (MPI_Offset*) malloc(nproc * sizeof(MPI_Offset));
    MPI_Allgather(&local_xsz, 1, MPI_OFFSET, all_collection_sizes, 1, MPI_OFFSET, MPI_COMM_WORLD);

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
    
    char* all_collections_buffer = (char*) malloc(total_recv_size);
    int* recvcounts =  (int*)malloc(nproc * sizeof(int));
    for (int i = 0; i < nproc; ++i) {
        recvcounts[i] = (int)all_collection_sizes[i];
    }
    
    MPI_Allgatherv(send_buffer, local_xsz, MPI_BYTE, all_collections_buffer, recvcounts, recv_displs, MPI_BYTE, MPI_COMM_WORLD);
    h5_grouparray **all_recv_meta = (h5_grouparray **)malloc(nproc * sizeof(h5_grouparray*));
    free_grouparray(local_meta);
    deserialize_all_grouparray(all_recv_meta, all_collections_buffer, recv_displs, recvcounts, nproc);
    // h5_grouparray * new_meta =  (h5_grouparray*)malloc(local_xsz);
    // deserialize_grouparray(new_meta, send_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    block_size = 2 * 1024 * 1024;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    H5Pset_coll_metadata_write(plist_id, true);
    H5Pset_meta_block_size(plist_id, block_size);
    outfile_id = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    create_all_metadata(all_recv_meta, nproc, outfile_id);

    // create_metadata(new_meta, outfile_id);
    io_time = MPI_Wtime() - end_time1;
    MPI_Barrier(MPI_COMM_WORLD);
    H5Pclose(plist_id);
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
    for (int i = 0; i < 5; i++){
        if (rank == 0) {
            printf("Max %s time: %f seconds\n", names[i], max_times[i]);
            printf("Min %s time: %f seconds\n", names[i], min_times[i]);
        }
    }
    // Free memory used by the group array
    free_all_grouparray(all_recv_meta, nproc);
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    free(recvcounts);
    MPI_Finalize();
    return 0;
}
