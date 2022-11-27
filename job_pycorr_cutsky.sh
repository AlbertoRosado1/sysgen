#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -L SCRATCH,project
#SBATCH -J pycorr
#SBATCH --mail-user=ar652820@ohio.edu
#SBATCH --mail-type=ALL
#SBATCH -o pycorr_zbin2_%J.out
#SBATCH -e pycorr_zbin2_%J.err
#SBATCH -t 08:00:00

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

srun -N 1 -n 1 python /global/homes/a/arosado/sysgen/pycorr_cutsky_all.py --zbin 2 --mockdir /global/cscratch1/sd/arosado/test_sysnet/null_mocks/