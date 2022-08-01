#!/bin/bash
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 5:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH,project
#SBATCH -J calculate_cl

# bash run_cl_random.sh

export PYTHONPATH=$HOME/sysnetdev/:$HOME/LSSutils/:$PYTHONPATH

module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
module use /global/common/software/m3169/cori/modulefiles
module load openmpi
module load python
conda activate sysnet

cd ${HOME}/sysgen/

do_cl=true

cl_random=${HOME}/sysgen/cl_random.py

if [ "${do_cl}" = true ]
then
    srun -n 1 python $cl_random
fi