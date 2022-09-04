import os
import numpy as np
import healpy as hp
from glob import glob
import fitsio as ft
import torch

loss_tol = -58.5

windows_l = []
windows_all = dict()
windows = dict()
for region in ['bmzls', 'ndecals', 'sdecals']:
    windows_all[region] = glob(f'/global/cscratch1/sd/arosado/test_sysnet/{region}_256/windows/window_model_*fits')
    for i, w_fn in enumerate(windows_all[region]):
        model_i = w_fn.split('_')[4]
        randseed = w_fn.split('_')[5] 
        snapshot_i = w_fn.split('_')[7].split('.')[0]
        snapshot = torch.load(f'/global/cscratch1/sd/arosado/test_sysnet/{region}_256/model_{model_i}_{randseed}/snapshot_{snapshot_i}.pth.tar')
        if snapshot['best_val_loss'] < loss_tol:
            windows_l.append(w_fn) 
    windows[region] = windows_l
    windows_l = []
    print(region)


nwindows = min([len(windows[r]) for r in ['bmzls', 'ndecals', 'sdecals']])
nside = 256

for i in range(nwindows):

    count_i = np.zeros(12*nside*nside)
    wind_i = np.zeros(12*nside*nside)
    
    for region in windows:
        d_ = ft.read(windows[region][i])
        wind_i[d_['hpix']] += d_['weight'] 
        count_i[d_['hpix']] += 1.0

    
    output_path = f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean/nnwindow_{i}.hp{nside}.fits'
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
         os.makedirs(output_dir)
    
    is_good = count_i > 0.0
    wind_i[is_good] = wind_i[is_good] / count_i[is_good]
    wind_i[~is_good] = hp.UNSEEN
    print(f'wrote {output_path}')
    hp.write_map(output_path, wind_i, dtype=np.float64, fits_IDL=False, overwrite=True)    
