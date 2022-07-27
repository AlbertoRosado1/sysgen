#!/bin/bash
# bash run_nnmean.sh

export PYTHONPATH=$HOME/sysnetdev/:$HOME/LSSutils/:$PYTHONPATH

module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
module use /global/common/software/m3169/cori/modulefiles
module load openmpi
module load python
conda activate sysnet

cd ${HOME}/sysgen/

do_nnmean=true

nnmean=${HOME}/sysgen/mean_nn_windows.py

if [ "${do_nnmean}" = true ]
then
    srun -n 1 python $nnmean
fi