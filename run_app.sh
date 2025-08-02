export PNETCDF_HINTS="nc_hash_size_dim=1048576;nc_hash_size_var=1048576"
mpiexec -n 4 ./app_baseline_test_all data/dataset_568k_metadata ./data/dataset_568k_metadata_app_baseline_test_all.h5
# mpiexec -n 4 ./app_baseline_test_all data/dataset_5m_metadata ./data/dataset_5m_metadata_app_baseline_test_all.nc
# read the created file
# mpiexec -n 4 ./app_baseline_read_test_all /files2/scratch/yll6162/parallel_metadata/app_baseline_test_all.nc