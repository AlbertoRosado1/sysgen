#!/bin/bash
cd $HOME

# for creating sysnet envirnoment on NERSC
module swap PrgEnv-${PE_ENV,,} PrgEnv-gnu
module use /global/common/software/m3169/cori/modulefiles
module load openmpi
module load python
conda create -n sysnet python=3.8 -y
conda activate sysnet

# custom build of mpi4py that is built on top of the OpenMPI loaded by the module
MPICC="mpicc -shared" pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

conda install scikit-learn jupyter ipykernel ipython matplotlib -y
conda install pytorch torchvision -c pytorch -y # cpu enabled 
conda install -c conda-forge fitsio healpy absl-py pytables pyyaml pandas -y # added pandas

#git clone https://github.com/mehdirezaie/sysnetdev.git # I already have cloned this
export PYTHONPATH=$HOME/sysnetdev/:$PYTHONPATH

python -m ipykernel install --user --name=sysnet --display-name "python (sysnet)"