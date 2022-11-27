import os
import utils as ut
import healpy as hp
import numpy as np
import fitsio as ft

norm_method = 'mean'
downsampling= 'frac'
#ph = 0
nside = 256
phs = [i for i in range(25)]
selection = dict(main=0, nz=0, Y5=1, sv3=0)

scratch = os.getenv('CSCRATCH')
out_dir = os.path.join(scratch,'test_sysnet','null_mocks')

# get NN-weights
nn_pred = hp.read_map(f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean_test1/nn-weights.hp256.fits')

for ph in phs:
    # Get null mock, we downsample it same way as contaminated mock but we shuffle the selection function, 
    # this way we have null mock with same dndz as cont mocks for fair comparison.
    mock_null = ut.get_mock(contaminated=True, selection_fn=nn_pred, nside=nside, norm_method=norm_method, downsampling=downsampling, seed=ph,
                            do_shuffle=True, tracer='LRG', ph=ph, **selection)

    fn = os.path.join(out_dir, f'cutsky_LRG_z0.800_AbacusSummit_base_c000_ph{ph:003d}.fits')
    ft.write(fn, mock_null)
    print(f"created null mock {fn}")