#!/bin/bash

#SBATCH -A m844
#SBATCH -t 00:5:00
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH --qos=debug

#SBATCH -o qout.256.%j # std::out 输出到这个文件 
#SBATCH -e qout.256.%j 

#SBATCH --mail-type=end,fail
#SBATCH --mail-user=youjia@northwestern.edu

if test "x$SLURM_NTASKS_PER_NODE" = x ; then #number of cores per node
   SLURM_NTASKS_PER_NODE=4 # 256 maximum
fi

NUM_NODES=$SLURM_JOB_NUM_NODES

NP=$((NUM_NODES * SLURM_NTASKS_PER_NODE))

ulimit -c unlimited

# EXE=/global/homes/y/yll6162/parallel_metadata/lib_baseline_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/save_input_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/app_baseline_read_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/app_baseline_test_all
EXE=/global/homes/y/yll6162/parallel_metadata/new_format_test_all
SB_EXE=/tmp/${USER}_test #tmp all users shared
sbcast -v ${EXE} ${SB_EXE} #executable copy to this position
# export LD_LIBRARY_PATH="/global/homes/y/yll6162/pnetcdf/pnetcdf-install/lib:$LD_LIBRARY_PATH"
# export PNETCDF_HINTS="nc_hash_size_dim=4096;nc_hash_size_var=4096"
# export PNETCDF_HINTS="nc_hash_size_dim=16777216;nc_hash_size_var=8388608"
srun -n $NP -c $((256/$SLURM_NTASKS_PER_NODE)) --cpu_bind=cores ${SB_EXE} 
#srun python ...