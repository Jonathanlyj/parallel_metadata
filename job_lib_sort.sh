#!/bin/bash

#SBATCH -A m844
#SBATCH -t 00:5:00
#SBATCH -N 4
#SBATCH -C cpu
#SBATCH --qos=debug

#SBATCH -o qout_log/lib_sort/qout_4n_64.256.%j # std::out is saved here
#SBATCH -e qout_log/lib_sort/qout_4n_64.256.%j 

#SBATCH --mail-type=end,fail
#SBATCH --mail-user=youjia@northwestern.edu

if test "x$SLURM_NTASKS_PER_NODE" = x ; then #number of cores per node
   SLURM_NTASKS_PER_NODE=64 # 256 max
fi

NUM_NODES=$SLURM_JOB_NUM_NODES

NP=$((NUM_NODES * SLURM_NTASKS_PER_NODE))

ulimit -c unlimited

EXE=/global/homes/y/yll6162/parallel_metadata/lib_baseline_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/save_input_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/app_baseline_read_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/app_baseline_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/new_format_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/mpi_io_test
# EXE=/global/homes/y/yll6162/parallel_metadata/free_timer
# EXE=/global/homes/y/yll6162/parallel_metadata/interleave_free_timer
# EXE=/global/homes/y/yll6162/parallel_metadata/sequential_free_timer
SB_EXE=/tmp/${USER}_test #tmp all users shared
sbcast -v ${EXE} ${SB_EXE} #executable copy to this position
# export LD_LIBRARY_PATH="/global/homes/y/yll6162/pnetcdf/pnetcdf-install/lib:$LD_LIBRARY_PATH"
# export PNETCDF_HINTS="nc_hash_size_dim=16777216;nc_hash_size_var=16777216"
export PNETCDF_HINTS="nc_hash_size_dim=16384;nc_hash_size_var=16384"

srun -n $NP -c $((256/$SLURM_NTASKS_PER_NODE)) --cpu_bind=cores ${SB_EXE} 
# srun -n $NP -c $((256/$SLURM_NTASKS_PER_NODE)) --cpu_bind=cores ${SB_EXE} -s 98 -c
#srun python ...