# parallel_metadata

## Run baseline method example on ECP
```
setenv PNETCDF_DIR /path/to/pnetcdf 

make

setenv LD_LIBRARY_PATH /path/to/pnetcdfPnetCDF-install/

e.g.
setenv LD_LIBRARY_PATH /files2/scratch/PnetCDF-meta/PnetCDF-install/lib

mpiexec -n {nproc} ./baseline_ex1


~/pnetcdf/PnetCDF-install/bin/ncmpidump lib_level_baseline.nc

```