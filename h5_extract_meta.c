#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include <inttypes.h>


#define SRC_FILE "nue_slice_panoptic_hdf_merged.h5"
// #define SRC_FILE "/files2/scratch/yll6162/parallel_metadata/script/test/nue_slice_graphs.0001.h5"
// #define SRC_FILE "h5_example.h5"
#define OUT_FILE "nue_slice_panoptic_hdf_merged_meta.h5"
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

        group_id = H5Gcreate2(file_id, group->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Iterate over each dataset in the group
        for (int j = 0; j < group->ndst; j++) {
            h5_dataset* dataset = group->datasets[j];

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

            dataset_id = H5Dcreate(group_id, dataset->name, datatype_id, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);


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
        H5Gclose(group_id);
    }
}



int main(int argc, char *argv[]) {

    hid_t outfile_id, plist_id ;
    h5_grouparray* local_meta;
    size_t local_xsz, global_xsz;

    local_meta = (h5_grouparray*)malloc(sizeof(h5_grouparray));
    
    read_metadata(local_meta);
    outfile_id = H5Fcreate(OUT_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    create_metadata(local_meta, outfile_id);

    H5Fclose(outfile_id);
    return 0;
}
