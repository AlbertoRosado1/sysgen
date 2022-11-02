#!/bin/bash
# bash run_nnsamp.sh

export PYTHONPATH=$HOME/sysnetdev/:$HOME/LSSutils/:$PYTHONPATH

module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
module use /global/common/software/m3169/cori/modulefiles
module load openmpi
module load python
conda activate sysnet

cd ${HOME}/sysgen/

do_nnsamp=true     # 3h x 10tpn

target=lrg
#region=$1   # options are bmzls, ndecals, sdecals
nside=256
version=test

nnsamp=${HOME}/sysgen/pull_sysnet_snapshot_mpidr9_mocks.py

if [ "${do_nnsamp}" = true ]
then
    srun -n 1 python $nnsamp #-n 1 python $nnsamp $region $version
fi