#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

#define FILE_NAME "dup.h5"
#define DATASET_NAME "my_dataset"

int main() {
    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[2] = {3, 3};
    
    // Create a new HDF5 file
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating file.\n");
        exit(1);
    }
    
    // Create the first dataset
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate2(file_id, DATASET_NAME, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset.\n");
        exit(1);
    }
    printf("Dataset '%s' created successfully.\n", DATASET_NAME);
    
    // Try to create a second dataset with the same name
    dataset_id = H5Dcreate2(file_id, DATASET_NAME, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // Close and release resources
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    
    return 0;
}
