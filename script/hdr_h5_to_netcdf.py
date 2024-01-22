import h5py
import pncpy
import os
zero_dim_ct = 0
dim_ct = 0

def copy_attributes(source, destination):
    for attr_name, attr_value in source.attrs.items():
        destination.put_att(attr_name, attr_value)


def hdf5_to_netcdf(hdf5_file, netcdf_file):
    # Open HDF5 file
    global zero_dim_ct
    global dim_ct
    with h5py.File(hdf5_file, 'r') as h5file:
        # Open NetCDF file for writing
        with pncpy.File(netcdf_file, 'w', format='64BIT_DATA') as ncfile:
            # TODO: Function to copy attributes

            # Iterate through HDF5 groups and datasets
            def convert_hdf5_to_netcdf(obj, parent_name=''):
                global zero_dim_ct
                global dim_ct
                for key, item in obj.items():
                    # Create NetCDF dimensions based on HDF5 dataset shape
                    if isinstance(item, h5py.Dataset):
                        item_name = parent_name + key
                        dims = []
                        for i in range(len(item.shape)):
                            dim_ct += 1
                            dim_name = f'{item_name}_dim_{i}'
                            if item.shape[i] != 0:
                                dim = ncfile.def_dim(dim_name, item.shape[i])
                                dims.append(dim)
                            else:
                                zero_dim_ct += 1
                        # Create NetCDF variables and copy data
                        # print(parent_name + key)
                        ncvar = ncfile.def_var(item_name, item.dtype, dimensions=tuple(dims))
                        # ncfile.enddef()
                        # ncvar[:] = item[:]
                        # Copy attributes from HDF5 to NetCDF variable
                        copy_attributes(item, ncvar)
                    elif isinstance(item, h5py.Group):
                        convert_hdf5_to_netcdf(item, f"{parent_name}_grp_{key}_")     
            convert_hdf5_to_netcdf(h5file)
    print(f"total dim: {dim_ct}, zero dim: {zero_dim_ct}")

if __name__ == "__main__":
    hdf5_file_path = "/files2/scratch/FNAL/uboone/pynuml_output/nue_slice_panoptic_hdf/nue_slice_graphs.0012.h5"
    nc_file_name = os.path.splitext(hdf5_file_path.split('/')[-1])[0] + ".nc"
    netcdf_file_path = nc_file_name
    hdf5_to_netcdf(hdf5_file_path, netcdf_file_path)
