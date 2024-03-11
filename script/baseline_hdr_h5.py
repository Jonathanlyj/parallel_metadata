import h5py
import numpy as np
import pickle
import mpi4py.MPI as MPI
import os
num_vars = 0
run_time_2 = 0
def create_example_hdf5_file(file_path):
    with h5py.File(file_path, 'w', libver='latest', driver='mpio', comm=MPI.COMM_WORLD) as file:
        # Dataset with Datatype, Dataspace, and Fill Value
        for i in range(5):
            group = file.create_group(f'group{i}')

            # Dataset with Datatype, Dataspace, and Fill Value inside the group
            dataset_group = group.create_dataset(f'dataset_inside_group{i}', compression='gzip', dtype='float64', shape=(10, 10), fillvalue=-1.0)
            # dataset_group = group.create_dataset(f'dataset_inside_group{i}', dtype='int32', shape=(10, 10), fillvalue=-1.0)

            # Adding Compression Information to the dataset inside the group
            # dataset_group.attrs.create('compression', 'gzip', dtype='S3')
            # dataset_group.attrs.create('compression_level', 9)

            # Adding Creation and Modification Time to the file
            file.attrs.create('creation_time', '2022-01-01 12:00:00')
            file.attrs.modify('modification_time', '2022-01-02 14:30:00')

def test_consistency_check(file_path, rank):
    with h5py.File(file_path, 'w', libver='latest', driver='mpio', comm=MPI.COMM_WORLD) as file:
        # Dataset with Datatype, Dataspace, and Fill Value
        for i in range(5):
            group = file.create_group(f'group{i}')

            # Dataset with Datatype, Dataspace, and Fill Value inside the group
            dataset_group = group.create_dataset(f'dataset_inside_group{i}_{rank}', compression='gzip', dtype='float64', shape=(10, 10), fillvalue=-1.0)
            # dataset_group = group.create_dataset(f'dataset_inside_group{i}', dtype='int32', shape=(10, 10), fillvalue=-1.0)

            # Adding Compression Information to the dataset inside the group
            # dataset_group.attrs.create('compression', 'gzip', dtype='S3')
            # dataset_group.attrs.create('compression_level', 9)

            # Adding Creation and Modification Time to the file
            file.attrs.create('creation_time', '2022-01-01 12:00:00')
            file.attrs.modify('modification_time', '2022-01-02 14:30:00')



# Traverse groups and datasets to display attributes



def store_metadata(instance, metadata):
    global num_vars
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
                num_vars += 1
                name_key = f'H5D_{name}'
            metadata[name_key] = {}
            store_metadata(item, metadata[name_key])

    elif isinstance(instance, h5py.Dataset):
        metadata['DTYPE'] = instance.dtype
        metadata['DSPACE'] = instance.shape

    return

def create_metadata(metadata, instance=None):
    global run_time_2

    for key, info in metadata.items():
        if key.startswith('H5A_'):
            # Attribute
            attr_name = key[len('H5A_'):]
            start_time_2 = MPI.Wtime()
            instance.attrs.create(attr_name, info['DAT'])
            end_time_2 = MPI.Wtime()
            run_time_2 += end_time_2 - start_time_2

        elif key.startswith('H5G_') or key.startswith('H5D_'):
            obj_name =  key[4:]
            if key.startswith('H5D_'):
                # Dataset
                dtype = info['DTYPE']
                dataspace = info['DSPACE']
                # dataset = instance.create_dataset(obj_name, shape = dataspace, dtype = dtype)
                start_time_2 = MPI.Wtime()
                plist = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
                plist.set_fill_time(h5py.h5d.FILL_TIME_NEVER)
                # if len(dataspace) > 1:
                #     chunk = tuple([1 for i in range(len(dataspace))])
                #     plist.set_chunk(chunk)
                h5py_dtype = h5py.h5t.py_create(dtype)
                spaceid = h5py.h5s.create_simple(dataspace)
                grpid = instance.id

                obj_name_b = bytes(obj_name,encoding="utf-8")
                datasetid = h5py.h5d.create(grpid, obj_name_b, h5py_dtype, spaceid, plist)
                dataset = h5py.Dataset(datasetid)
                create_metadata(info, instance=dataset)
                end_time_2 = MPI.Wtime()
                run_time_2 += end_time_2 - start_time_2
            else:
                # Group
                start_time_2 = MPI.Wtime()
                group = instance.create_group(obj_name)
                end_time_2 = MPI.Wtime()
                run_time_2 += end_time_2 - start_time_2
                create_metadata(info, instance=group)

def create_file_from_metadata(file_path, metadata, meta_block_size=None):
    with h5py.File(file_path, 'w', driver="mpio", libver='latest', comm=MPI.COMM_WORLD, meta_block_size=meta_block_size) as file:
        create_metadata(metadata, instance=file)


def store_file_metadata(file_path):
    metadata = {}
    with h5py.File(file_path, 'r') as file:
        # for key, value in file['group2']['dataset_inside_group2'].attrs.items():
        #     print(key, value)
        store_metadata(file, metadata)
    return metadata

def store_file_metadata_parallel(file_paths, rank, size):
    global num_vars
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
                        num_vars += 1
                        name_key = f'H5D_{name}'
                    metadata[name_key] = {}
                    store_metadata(item, metadata[name_key])
    return metadata


def exchange_metadata(metadata, rank, size):
    global num_vars
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
        print(f"rank 0 n_vars: {num_vars}")
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
    # app_file_dir = "/homes/yll6162/parallel_metadata/script/test"

    mb = 1024 * 1024
    meta_block_size = 16
    app_h5_files = []
    for filename in os.listdir(app_file_dir):
    # Check if the file has a ".h5" extension
        if filename.endswith(".h5"):
            # Construct the full file path
            file_path = os.path.join(app_file_dir, filename)
            app_h5_files.append(file_path)

    dummy_file_path = 'dummy_test.h5'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # test_consistency_check(dummy_file_path, rank)
    # ---------------------------------------------------- Read Metadata ----------------------------------------------------
    # create_example_hdf5_file(dummy_file_path)
    metadata = store_file_metadata_parallel(app_h5_files[:], rank, size)

    # ---------------------------------------------------- Exchange Metadata ------------------------------------------------
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
        

    # ---------------------------------------------------- Create Metadata Collectively--------------------------------------------------
    output_file_path = f"{app_file_dir.split('/')[-1]}_merged.h5"
    comm.Barrier()
    # start_time_2 = MPI.Wtime()
    if meta_block_size:
        create_file_from_metadata(output_file_path, all_metadata, meta_block_size = meta_block_size * mb)
    else:
        create_file_from_metadata(output_file_path, all_metadata)
    # end_time_2 = MPI.Wtime()
    # crt_time = end_time_2 - start_time_2
    min_time = comm.reduce(run_time_2, op=MPI.MIN, root=0)
    max_time = comm.reduce(run_time_2, op=MPI.MAX, root=0)
    if rank==0:
        print(f"Max create time: {max_time}")
        print(f"Min create time: {min_time}")
        print(f"Meta block size set: {meta_block_size} MB")

    


