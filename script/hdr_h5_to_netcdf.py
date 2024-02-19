import h5py
import pncpy
import os
import argparse
zero_dim_ct = 0
dim_ct = 0
var_ct = 0
def copy_attributes(source, destination):
    for attr_name, attr_value in source.attrs.items():
        destination.put_att(attr_name, attr_value)


def hdf5_to_netcdf(hdf5_file, netcdf_file, use_small):
    # Open HDF5 file
    global zero_dim_ct
    global dim_ct
    global var_ct
    with h5py.File(hdf5_file, 'r') as h5file:
        # Open NetCDF file for writing
        with pncpy.File(netcdf_file, 'w', format='64BIT_DATA') as ncfile:
            # TODO: Function to copy attributes

            # Iterate through HDF5 groups and datasets
            def convert_hdf5_to_netcdf(obj, parent_name=''):
                global zero_dim_ct
                global dim_ct
                global var_ct
                for key, item in obj.items():
                    # Create NetCDF dimensions based on HDF5 dataset shape
                    if isinstance(item, h5py.Dataset):
                        item_name = parent_name + key
                        var_ct += 1
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
                        if not use_small or hash(key) % 10 == 0:
                            convert_hdf5_to_netcdf(item, f"{parent_name}_grp_{key}_")

            convert_hdf5_to_netcdf(h5file)
            ncfile.enddef()
    print(f"total dim: {dim_ct}, total variable: {var_ct}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert HDF5 groups to NetCDF")
    parser.add_argument("--input_file", help="Path to the input HDF5 file",
                        default="/files2/scratch/FNAL/uboone/pynuml_output/nue_slice_panoptic_hdf/nue_slice_graphs.0001.h5")
    parser.add_argument("--output_file", help="Path to the output NetCDF file",
                        default="./nue_slice_graphs.0001.nc")
    parser.add_argument("--small", action="store_true", help="Use 1/10 th of groups")

    args = parser.parse_args()
    use_small = False
    if args.small:
        use_small = True
    hdf5_to_netcdf(args.input_file, args.output_file, use_small)
