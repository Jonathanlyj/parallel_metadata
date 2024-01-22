# parallel_metadata

## Run baseline method example on ECP
```
setenv PNETCDF_DIR /path/to/pnetcdf 

make

setenv LD_LIBRARY_PATH /path/to/pnetcdfPnetCDF-install/

mpiexec -n {nproc} ./baseline_ex1

```