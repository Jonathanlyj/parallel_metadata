#include "hdf5.h"

#define FILE_NAME "h5_example.h5"
#define GROUP_NAMES {"Group1", "Group2", "Group3", "Group4"}
#define DATASET_NAME "Dataset1"

int main(void) {
    hid_t file, group, dataset, dataspace, dcpl;
    hsize_t dims[2] = {2, 3};
    herr_t             status;
    char *group_names[] = GROUP_NAMES;

    // Create a new HDF5 file
    file = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create groups under the root
    for (int i = 0; i < 4; i++) {
        group = H5Gcreate2(file, group_names[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create a dataspace
        dataspace = H5Screate_simple(2, dims, NULL);
        dcpl = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_alloc_time(dcpl, H5D_ALLOC_TIME_LATE);
        // Create a dataset
        dataset = H5Dcreate2(group, DATASET_NAME, H5T_IEEE_F32LE, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

        // Write different data to each dataset
        int data[2][3];
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 3; k++) {
                data[j][k] = (i + 1) * (j + 1) * (k + 1); // Different formula for each group
            }
        }
        // H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

        // Close resources for this group
        H5Pclose(dcpl);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Gclose(group);
    }

    // Close the file
    H5Fclose(file);

    return 0;
}
