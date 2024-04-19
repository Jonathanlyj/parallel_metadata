#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hdf5.h"

// #define SRC_FILE "nue_slice_panoptic_hdf_merged.h5"
#define SRC_FILE "h5_example.h5"
#define OUT_FILE "h5_baseline_test_all.h5"
double crt_start_time, total_crt_time=0;
// Define a structure to represent a dataset
typedef struct {
    char *name;
    int name_len;
    hid_t datatype;
    int ndims;
    hsize_t *dims;
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

// Function to create a group
void new_group(const char *name, int name_len, int ndst, h5_group* group) {
    group->name = strdup(name);
    group->name_len = name_len;
    group->ndst = ndst;
    group->datasets = (h5_dataset**)malloc(ndst * sizeof(h5_dataset*));
    return;
}


void new_dataset(const char *name, int name_len, hid_t datatype, int ndims, hsize_t* dims, h5_dataset* dataset){
    dataset->name = strdup(name);
    dataset->name_len = name_len;
    dataset->ndims = ndims;
    dataset->datatype = datatype;
    dataset->dims = malloc(ndims * sizeof(hsize_t));
    for (int i = 0; i < ndims; i++){
        dataset->dims[i] = dims[i];
    }
}

// Function to free memory used by a group
void free_dataset(h5_dataset *dataset) {
    if (dataset == NULL) {
        return;
    }
    free(dataset->name);
    free(dataset->dims);
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
    hid_t       file_id, group_id, dst_id, datatype, dataspace, memspace;
    herr_t      status;
    hsize_t ngrp, ndst;
    hsize_t *dims_out;
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
            H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);

            dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
            // Get dataspace and datatype
            dataspace = H5Dget_space(dst_id);
            datatype = H5Dget_type(dst_id);

            // Get dataset dimensions
            int ndims;
            ndims = H5Sget_simple_extent_ndims(dataspace);
            dims_out = (hsize_t *)malloc(ndims * sizeof(hsize_t));
            status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
            new_dataset(dst_name, strlen(dst_name),  datatype, ndims, dims_out, local_meta->groups[i]->datasets[j]);
            free(dims_out);
            H5Dclose(dst_id);
        }
        // Close the group
        H5Gclose(group_id);
    }
    // Close/release resources
    H5Fclose(file_id);
}

void read_metadata_parallel(h5_grouparray* local_meta, int rank, int nproc){
    hid_t       file_id, group_id, dst_id, datatype, dataspace, memspace;
    herr_t      status;
    hsize_t ngrp, ndst;
    hsize_t *dims_out;
    int grps_per_proc, remainder, start, count;



   // Open the file
    file_id = H5Fopen(SRC_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

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
            local_meta->groups[i - start]->datasets[j] = malloc(sizeof(h5_dataset));
            H5Gget_objname_by_idx(group_id, (hsize_t)j, dst_name, 256);

            dst_id = H5Dopen2(group_id, dst_name, H5P_DEFAULT); 
            // Get dataspace and datatype
            dataspace = H5Dget_space(dst_id);
            datatype = H5Dget_type(dst_id);

            // Get dataset dimensions
            int ndims;
            ndims = H5Sget_simple_extent_ndims(dataspace);
            dims_out = (hsize_t *)malloc(ndims * sizeof(hsize_t));
            status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
            new_dataset(dst_name, strlen(dst_name),  datatype, ndims, dims_out, local_meta->groups[i - start]->datasets[j]);
            free(dims_out);
            H5Dclose(dst_id);
        }
        // Close the group
        H5Gclose(group_id);
    }
    // Close/release resources
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
            total_size += sizeof(hid_t); // For datatype
            total_size += sizeof(int); // For ndims
            total_size += dataset->ndims * sizeof(hsize_t); // For dims
        }
    }
    return total_size;
}
void serialize_dataset(h5_dataset *dataset, char **ptr) {
    memcpy(*ptr, &(dataset->name_len), sizeof(int));
    *ptr += sizeof(int);
    memcpy(*ptr, dataset->name, dataset->name_len);
    *ptr += dataset->name_len;
    memcpy(*ptr, &(dataset->datatype), sizeof(hid_t));
    *ptr += sizeof(hid_t);
    memcpy(*ptr, &(dataset->ndims), sizeof(int));
    *ptr += sizeof(int);
    memcpy(*ptr, dataset->dims, dataset->ndims * sizeof(hsize_t));
    *ptr += dataset->ndims * sizeof(hsize_t);
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
    memcpy(&(dataset->datatype), *ptr, sizeof(hid_t));
    *ptr += sizeof(hid_t);
    memcpy(&(dataset->ndims), *ptr, sizeof(int));
    *ptr += sizeof(int);
    dataset->dims = (hsize_t *)malloc(dataset->ndims * sizeof(hsize_t));
    memcpy(dataset->dims, *ptr, dataset->ndims * sizeof(hsize_t));
    *ptr += dataset->ndims * sizeof(hsize_t);
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

static int deserialize_all_grouparray(h5_grouparray **all_recv_hdr, char* all_collections_buffer, int* recv_displs, int* recvcounts, int nproc){
    for (int i=0; i< nproc; i++){
        all_recv_hdr[i]= (h5_grouparray *)malloc(sizeof(h5_grouparray));
        deserialize_grouparray(all_recv_hdr[i], all_collections_buffer + recv_displs[i]);
    }
    return 0;
}

static int free_all_grouparray(h5_grouparray **all_recv_hdr, int nproc){
    if (all_recv_hdr != NULL){
        for (int i=0; i< nproc; i++) free_grouparray(all_recv_hdr[i]);
        free(all_recv_hdr);
    }
    return 0;
}

/* ---------------------------------- Write Metadata ----------------------------------------*/
void create_metadata(h5_grouparray* local_meta, hid_t file_id) {
    hid_t group_id, dataset_id, dataspace, datatype;
    herr_t status;

    // Iterate over each group in the local_meta
    for (int i = 0; i < local_meta->ngrps; i++) {
        h5_group* group = local_meta->groups[i];

        // Create a group in the HDF5 file
        crt_start_time = MPI_Wtime();
        group_id = H5Gcreate2(file_id, group->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        total_crt_time += MPI_Wtime() - crt_start_time;

        // Iterate over each dataset in the group
        for (int j = 0; j < group->ndst; j++) {
            h5_dataset* dataset = group->datasets[j];

            // Create dataspace
            crt_start_time = MPI_Wtime();
            dataspace = H5Screate_simple(dataset->ndims, dataset->dims, NULL);
            total_crt_time += MPI_Wtime() - crt_start_time;
            // Create datatype
            datatype = dataset->datatype;
            // Create dataset
            crt_start_time = MPI_Wtime();
            dataset_id = H5Dcreate2(group_id, dataset->name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            total_crt_time += MPI_Wtime() - crt_start_time;

            // Close resources
            H5Sclose(dataspace);
            H5Dclose(dataset_id);
        }

        // Close the group
        H5Gclose(group_id);
    }
}

int create_all_metadata(h5_grouparray **all_recv_hdr, int nproc, hid_t file_id){
    for (int i=0; i< nproc; i++){
        h5_grouparray *hdr_data = all_recv_hdr[i];
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

    // status = serialize_meta(&local_meta, send_buffer);
    // // Access group and dataset information
    // for (int i = 0; i < local_meta->ngrps; i++) {
    //     printf("Group name: %s, Length: %d\n", local_meta->groups[i]->name, local_meta->groups[i]->name_len);
    //     printf("Number of datasets: %d\n", local_meta->groups[i]->ndst);
    //     printf("First dataset name: %s, Length: %d\n", local_meta->groups[i]->datasets[0]->name, local_meta->groups[i]->datasets[0]->name_len);
    // }
       // printf("rank %d, buffer size: %lld \n", rank, local_hdr.xsz);


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
    h5_grouparray **all_recv_hdr = (h5_grouparray **)malloc(nproc * sizeof(h5_grouparray*));
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    outfile_id = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    deserialize_all_grouparray(all_recv_hdr, all_collections_buffer, recv_displs, recvcounts, nproc);
    MPI_Barrier(MPI_COMM_WORLD);
    end_time1 = MPI_Wtime();
    create_all_metadata(all_recv_hdr, nproc, outfile_id);

    MPI_Barrier(MPI_COMM_WORLD);
    end_time3 = MPI_Wtime();
    H5Pclose(plist_id);
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
    free_all_grouparray(all_recv_hdr, nproc);
    free_grouparray(local_meta);
    free(send_buffer);
    free(all_collections_buffer);
    free(all_collection_sizes);
    free(recv_displs);
    free(recvcounts);
    MPI_Finalize();
    return 0;
}