import h5py
import numpy as np
import pickle
import mpi4py.MPI as MPI
import pncpy as pnc
import os
import sys

ncopy = 10

zero_dim_ct = 0
dim_ct = 0
var_ct = 0
start_time_1 = 0

# run_time_2 = 0
def test_consistency_check(file_path, rank, comm):

    file = pnc.File(file_path, 'w', comm=comm,  format='64BIT_DATA')
        # Dataset with Datatype, Dataspace, and Fill Value
    for i in range(1):
        # Dataset with Datatype, Dataspace, and Fill Value inside the group
        if rank < 3:
            dim = file.def_dim(f'dataset_{rank}_dim_{i}', 10)
            print(f"rank{rank}: def dim completed")
            sys.stdout.flush()
            var = file.def_var(f'dataset_{rank}_{i}', np.dtype(int), dimensions = [dim])
            print(f"rank{rank}: def var completed")
            sys.stdout.flush()
        # dataset_group = group.create_dataset(f'dataset_inside_group{i}', dtype='int32', shape=(10, 10), fillvalue=-1.0)
    print(f"rank{rank}: before enddef")
    sys.stdout.flush()
    # flush
    # command
    file.enddef()
    # file.close()
    print(f"rank{rank}: after enddef")
    sys.stdout.flush()
    file.close()
    comm.Barrier()
    print(f"rank{rank}: after file close")
    # sys.stdout.flush()
def create_metadata_nc(metadata, instance, parent_name=''):
    global zero_dim_ct
    global dim_ct
    global var_ct
    global run_time_2
    for key, info in metadata.items():
        # Create NetCDF dimensions based on HDF5 dataset shape
        if key.startswith('H5A_'):
            # Attribute
            attr_name = key[len('H5A_'):]
            # start_time_2 = MPI.Wtime()
            instance.put_att(attr_name, info['DAT'])
            # end_time_2 = MPI.Wtime()
            # run_time_2 += end_time_2 - start_time_2
        elif key.startswith('H5G_') or key.startswith('H5D_'):
            
            obj_name =  key[4:]
            item_name = parent_name + obj_name
            
            if key.startswith('H5D_'):
                dtype = info['DTYPE']
                dataspace = info['DSPACE']
                dims = []
            
                for i in range(len(dataspace)):
                    dim_ct += 1
                    dim_name = f'{item_name}_dim_{i}'
                    if dataspace[i] != 0:
                        # start_time_2 = MPI.Wtime()
                        dim = instance.def_dim(dim_name, dataspace[i])
                        # end_time_2 =  MPI.Wtime()
                        # run_time_2 += end_time_2 - start_time_2
                        dims.append(dim)
                    else:
                        zero_dim_ct += 1
                # Create NetCDF variables and copy data
                # start_time_2 = MPI.Wtime()
                var_ct += 1
                ncvar = instance.def_var(item_name, dtype, dimensions=tuple(dims))
                # end_time_2 =  MPI.Wtime()
                # run_time_2 += end_time_2 - start_time_2
                create_metadata_nc(info, ncvar,item_name)
            else:
            #group
                create_metadata_nc(info, instance, f"{parent_name}_grp_{obj_name}_")






# Traverse groups and datasets to display attributes



def store_metadata(instance, metadata):
    # global var_ct
    for key, value in instance.attrs.items():
        att_key = f'H5A_{key}'
        metadata[att_key] = {}
        metadata[att_key]['DAT'] = value
        if not isinstance(instance.attrs[key], str):
            metadata[att_key]['DTYPE'] =  value.dtype
            metadata[att_key]['DSPACE'] = value.shape

    if isinstance(instance, h5py.Group):
        for name, item in instance.items():
            if isinstance(item, h5py.Group):
                name_key = f'H5G_{name}'
            elif isinstance(item, h5py.Dataset):
                # var_ct += 1
                name_key = f'H5D_{name}'
            metadata[name_key] = {}
            store_metadata(item, metadata[name_key])

    elif isinstance(instance, h5py.Dataset):
        metadata['DTYPE'] = instance.dtype
        metadata['DSPACE'] = instance.shape

    return


def create_ncfile_from_metadata(file_path, metadata, ncopy):
    global start_time_2
    global run_time_2

    with pnc.File(file_path, 'w', comm=MPI.COMM_WORLD,  format='64BIT_DATA') as file:
        start_time_2 = MPI.Wtime()
        for i in range(ncopy):
            create_metadata_nc(metadata, instance=file, parent_name=f'copy_{i}_')
        # start_time_2 = MPI.Wtime()
        file.enddef()
        # end_time_2 = MPI.Wtime()
        # run_time_2 += end_time_2 - start_time_2



def store_file_metadata(file_path):
    metadata = {}
    with h5py.File(file_path, 'r') as file:
        # for key, value in file['group2']['dataset_inside_group2'].attrs.items():
        #     print(key, value)
        store_metadata(file, metadata)
    return metadata

def store_file_metadata_parallel(file_paths, rank, size):
    global var_ct
    metadata = {}
    # Only rank 0 stores file-level attributes and file-level datasets
    for file_idx, file_path in enumerate(file_paths):
        filename = file_path.split('/')[-1].rstrip('.h5')
        with h5py.File(file_path, 'r', libver='latest', driver = 'mpio', comm=MPI.COMM_WORLD) as file:
            for key, value in file.attrs.items():
                att_key = f'H5A_{filename}_{key}'
                if rank == file_idx % size:
                    metadata[att_key] = {}
                    metadata[att_key]['DAT'] = value
                    if not isinstance(value, str):
                        metadata[att_key]['DTYPE'] =  value.dtype
                        metadata[att_key]['DSPACE'] = value.shape
        
            # Partition other metadata by root-level groups/datasets
            for idx, (name, item) in enumerate(file.items()):
                if rank == idx % size:
                    if isinstance(item, h5py.Group):
                        name_key = f'H5G_{name}'
                    elif isinstance(item, h5py.Dataset):
                        # var_ct += 1
                        name_key = f'H5D_{name}'
                    metadata[name_key] = {}
                    store_metadata(item, metadata[name_key])
    return metadata


def exchange_metadata(metadata, rank, size):
    global var_ct

    comm = MPI.COMM_WORLD
    # Serialize
    buffer = pickle.dumps(metadata)
    send_count = len(buffer)

    send_counts = comm.allgather(send_count)
    # Determine the displacement for each process
    send_displacements = np.cumsum(send_counts) - send_counts
    # Gather data from all processes
    all_metadata_buff = np.empty(sum(send_counts), dtype='S1')
    if rank == 0: 
        print(f"Total size: {sum(send_counts)/(1024 * 1024)} MB")
        print(f"Max size: {max(send_counts)/(1024 * 1024)} MB")
        print(f"Min size: {min(send_counts)/(1024 * 1024)} MB")
        print(f"rank 0 n_vars: {var_ct}")
    comm.Allgatherv(buffer, [all_metadata_buff, send_counts, send_displacements, MPI.CHAR])
    # Deserialize
    all_metadata = {}
    for i in range(size):
        meta = pickle.loads(all_metadata_buff[send_displacements[i]:send_displacements[i] + send_counts[i]])
        all_metadata.update(meta)
    return all_metadata





if __name__ == "__main__":
    # app_file_path = "/files2/scratch/FNAL/uboone/pynuml_output/nue_slice_panoptic_hdf/nue_slice_graphs.0001.h5"
    app_file_dir = "/files2/scratch/FNAL/uboone/pynuml_output/nue_slice_panoptic_hdf"
    # app_file_dir = "./test"
    
    app_h5_files = []
    for filename in os.listdir(app_file_dir):
    # Check if the file has a ".h5" extension
        if filename.endswith(".h5"):
            # Construct the full file path
            file_path = os.path.join(app_file_dir, filename)
            app_h5_files.append(file_path)

    dummy_file_path = 'dummy_test.nc'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # test_consistency_check(dummy_file_path, rank, comm)
    # ---------------------------------------------------- Read Metadata ----------------------------------------------------
    # create_example_hdf5_file(dummy_file_path)
    metadata = store_file_metadata_parallel(app_h5_files[:], rank, size)
    

    # # ---------------------------------------------------- Exchange Metadata ------------------------------------------------
    comm.Barrier()
    start_time_1 = MPI.Wtime()
    all_metadata = exchange_metadata(metadata, rank, size)
    end_time_1 = MPI.Wtime()
    comm_time = end_time_1 - start_time_1
    min_time = comm.reduce(comm_time, op=MPI.MIN, root=0)
    max_time = comm.reduce(comm_time, op=MPI.MAX, root=0)
    if rank==0:
        print(f"Max comm time: {max_time}")
        print(f"Min comm time: {min_time}")

    # # ---------------------------------------------------- Create Metadata Collectively--------------------------------------------------

    output_file_path = f"{app_file_dir.split('/')[-1]}_merged_{ncopy}_copy.nc"

    create_ncfile_from_metadata(output_file_path, all_metadata, ncopy)
    end_time_2 = MPI.Wtime()
    crt_time = end_time_2 - start_time_2
    min_time = comm.reduce(crt_time, op=MPI.MIN, root=0)
    max_time = comm.reduce(crt_time, op=MPI.MAX, root=0)
    if rank==0:
        print(f"Max create time: {max_time}")
        print(f"Min create time: {min_time}")
        print(f"total dim: {dim_ct}, total variable: {var_ct}")


