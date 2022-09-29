import os
import numpy as np
import healpy as hp
from glob import glob
import fitsio as ft
import torch

#loss_tol = [-58.5, -63.1, -61.2]
t = '_test1'
bins_h = 100
regions = ['bmzls', 'ndecals', 'sdecals']
snapshots = dict()
loss_l = []
loss =  dict()
for region in regions:
    snapshots[region] = glob(f'/global/cscratch1/sd/arosado/test_sysnet/{region}_256/model_*_*/snapshot_*.pth.tar')
    for i, sn_fn in enumerate(snapshots[region]):
        snapshot = torch.load(sn_fn)
        loss_l.append(snapshot['best_val_loss'])
    loss[region] = loss_l
    loss_l = []
loss_tol = [((np.histogram(loss[r], bins=bins_h)[1][:-1] + np.histogram(loss[r], bins=bins_h)[1][1:])/2.).min() for r in regions]
print(f'loss_tol {regions}: {loss_tol}')

windows_l = []
windows_all = dict()
windows = dict()
for loss_i, region in enumerate(regions):
    windows_all[region] = glob(f'/global/cscratch1/sd/arosado/test_sysnet/{region}_256/windows/window_model_*fits')
    for i, w_fn in enumerate(windows_all[region]):
        model_i = w_fn.split('_')[4]
        randseed = w_fn.split('_')[5] 
        snapshot_i = w_fn.split('_')[7].split('.')[0]
        snapshot = torch.load(f'/global/cscratch1/sd/arosado/test_sysnet/{region}_256/model_{model_i}_{randseed}/snapshot_{snapshot_i}.pth.tar')
        loss_l.append(snapshot['best_val_loss'])
        if snapshot['best_val_loss'] < loss_tol[loss_i]:
            windows_l.append(w_fn) 
    windows[region] = windows_l
    windows_l = []
    print(region)


nwindows = min([len(windows[r]) for r in regions])
print(f'nwindows: {nwindows}')

nside = 256

for i in range(nwindows):

    count_i = np.zeros(12*nside*nside)
    wind_i = np.zeros(12*nside*nside)
    
    for region in windows:
        d_ = ft.read(windows[region][i])
        wind_i[d_['hpix']] += d_['weight'] 
        count_i[d_['hpix']] += 1.0

    
    output_path = f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean{t}/nnwindow_{i}.hp{nside}.fits'
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
         os.makedirs(output_dir)
    
    is_good = count_i > 0.0
    wind_i[is_good] = wind_i[is_good] / count_i[is_good]
    wind_i[~is_good] = hp.UNSEEN
    print(f'wrote {output_path}')
    hp.write_map(output_path, wind_i, dtype=np.float64, fits_IDL=False, overwrite=True)    
