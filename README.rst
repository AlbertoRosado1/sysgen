sysgen
======

Repo for generating systematics on NERSC using Mehdi Rezaie's NN (SYSNet).

Step 1: Train the NN with real data. (repeat for each region)
--------

Create interactive job and inside ``run_sysnet.sh`` set ``do_LRfinder=true``, then run::

    salloc -N 1 -C haswell -t 04:00:00 --qos interactive -L SCRATCH,project -J sysgen
    bash run_sysnet.sh $region
    
.. note:: options for regions are: bmzls, ndecals, sdecals

Look for the directory that contains ``loss_vs_lr_0.png`` and use this plot to determine the learning rate for the NN. Set this learning rate (e.g. ``lr=0.005``) and set ``do_LRfinder=false`` and ``do_nnrun=true`` . Now run::

    bash run_sysnet.sh $region
    
Step 2 - part 1: Forward pass given the models obtained by training the NN. (repeat for each region)
--------

Create interactive job (first ``exit`` the previous one if still have time remianing)::

    salloc -N 1 -C haswell -t 04:00:00 --qos interactive -L SCRATCH,project -J sysgen
