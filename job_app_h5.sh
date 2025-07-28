#!/bin/bash

#SBATCH -A m844
#SBATCH -t 00:29:00
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH --qos=debug

#SBATCH -o qout_log/app_h5/qout_1n_64.256.%j # std::out is saved here
#SBATCH -e qout_log/app_h5/qout_1n_64.256.%j 

#SBATCH --mail-type=end,fail
#SBATCH --mail-user=youjia@northwestern.edu

if test "x$SLURM_NTASKS_PER_NODE" = x ; then #number of cores per node
   SLURM_NTASKS_PER_NODE=64 # 256 max
fi

NUM_NODES=$SLURM_JOB_NUM_NODES

NP=$((NUM_NODES * SLURM_NTASKS_PER_NODE))

ulimit -c unlimited

# EXE=/global/homes/y/yll6162/parallel_metadata/lib_baseline_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/save_input_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/app_baseline_read_test_all
EXE=/global/homes/y/yll6162/parallel_metadata/h5_baseline_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/new_format_test_all
# EXE=/global/homes/y/yll6162/parallel_metadata/mpi_io_test
# EXE=/global/homes/y/yll6162/parallel_metadata/free_timer
# EXE=/global/homes/y/yll6162/parallel_metadata/interleave_free_timer
# EXE=/global/homes/y/yll6162/parallel_metadata/sequential_free_timer
SB_EXE=/tmp/${USER}_test #tmp all users shared
sbcast -v ${EXE} ${SB_EXE} #executable copy to this position


srun -n $NP -c $((256/$SLURM_NTASKS_PER_NODE)) --cpu_bind=cores ${SB_EXE} 
# srun -n $NP -c $((256/$SLURM_NTASKS_PER_NODE)) --cpu_bind=cores ${SB_EXE} -s 98 -c
#srun python ...