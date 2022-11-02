#!/bin/bash
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 05:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH,project
#SBATCH -J create_selection_func_ens

run_sysnet=/global/homes/a/arosado/sysnetdev/scripts/app.py

export PYTHONPATH=$HOME/sysnetdev/:$HOME/LSSutils/:$PYTHONPATH

module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
module use /global/common/software/m3169/cori/modulefiles
module load openmpi
module load python
conda activate sysnet

# nn parameters
axes=(0 4 10) #  extinction, galdepth-z, and psfsize-r # (0 1 2 3 4 5 6 7 8 9 10 11 12) 
nchain=5
nepoch=200
nns=(4 20)
bsize=5000
lr=0.005
model=dnnp
loss=pnll
etamin=0.00001

do_LRfinder=false # for running the learning rate finder
do_nnrun=true

if [ "${do_LRfinder}" = true ]
then
       input_data=/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/cutsky_LRG_z0.800_AbacusSummit_base_c000_ph000.fits
       output_nn=/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/
       du -h $input_data
       srun -n 1 python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -fl
fi

if [ "${do_nnrun}" = true ]
then
       input_data=/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/cutsky_LRG_z0.800_AbacusSummit_base_c000_ph000.fits
       output_nn=/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/
       du -h $input_data
       srun -n 1 python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -lr $lr --eta_min $etamin -ne $nepoch -nc $nchain --snapshot_ensemble -k --no_eval
fi