#!/bin/bash
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 20:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH,project
#SBATCH -J create_selection_func_ens

run_sysnet=/global/homes/a/arosado/sysnetdev/scripts/appensemble.py

export PYTHONPATH=$HOME/sysnetdev/:$PYTHONPATH
conda activate sysnet

# regions
region=$1 # options are bmzls, ndecals, sdecals

# nn parameters
axes=(0 1 2 3 4 5 6 7 8 9 10 11 12) 
nchain=5
nepoch=200
nns=(4 20)
bsize=5000
lr=0.005
model=dnnp
loss=pnll
etamin=0.00001

do_LRfinder=true # for running the learning rate finder
do_nnrun=false

if [ "${do_LRfinder}" = true ]
then
       input_data=/global/cscratch1/sd/arosado/nlrg_features_${region}_256.fits
       output_nn=/global/cscratch1/sd/arosado/test_sysnet/${region}_256/
       du -h $input_data
       srun -n 1 python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -fl
fi

if [ "${do_nnrun}" = true ]
then
       input_data=/global/cscratch1/sd/arosado/nlrg_features_${region}_256.fits
       output_nn=/global/cscratch1/sd/arosado/test_sysnet/${region}_256/
       du -h $input_data
       srun -n 1 python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -lr $lr --eta_min $etamin -ne $nepoch -nc $nchain --snapshot_ensemble -k --no_eval
fi