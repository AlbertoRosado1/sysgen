sysgen
======

Repo for generating systematics on NERSC using Mehdi Rezaie's NN (SYSNet).

For running this code you will need to run ``create_sysnet_env.sh`` at least once (assuming you have not approppiately installed SYSNet before). Also you will need to clone https://github.com/mehdirezaie/LSSutils.

Step 1: Train the NN with real data. (repeat for each region)
--------

Create interactive job and inside ``run_sysnet.sh`` set ``do_LRfinder=true``, then run::

    salloc -N 1 -C haswell -t 04:00:00 --qos interactive -L SCRATCH,project -J sysgen
    bash run_sysnet.sh $region
    
.. note:: options for regions are: bmzls, ndecals, sdecals

Look for the directory that contains ``loss_vs_lr_0.png`` and use this plot to determine the learning rate for the NN. Set this learning rate (e.g. ``lr=0.005``) and set ``do_LRfinder=false`` and ``do_nnrun=true`` . Now run::

    bash run_sysnet.sh $region
    exit
    
Step 2 - part 1: Forward pass given the models obtained by training the NN. (repeat for each region)
--------

Create new interactive job (first ``exit`` the previous if there is time allocation remaining)::

    salloc -N 1 -C haswell -t 04:00:00 --qos interactive -L SCRATCH,project -J sysgen
    bash run_nnsamp.sh $region 

Step 2 - part 2: Combine all windows obtained from forward pass
--------

If the interactive job time allocation has expired create another and run::

    bash run_nncomb.sh

Then run ``run_nnmean.sh`` to obtain one selection function from the mean of the ensemble of selection functions::
    
    bash run_nnmean.sh
    exit
   
Step 3: Contaminate mocks (subsampling)
--------

Look at ``subsample_mock_nb.ipynb`` notebook.

