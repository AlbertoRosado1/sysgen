import os
import numpy as np
import healpy as hp
from glob import glob
import fitsio as ft
import torch

#loss_tol = [-58.5, -63.1, -61.2]
t = ''#'_test1'

windows = glob(f'/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/windows/window_model_*fits')

nwindows = len(windows)
print(f'nwindows: {nwindows}')

nside = 256

for i in range(nwindows):

    count_i = np.zeros(12*nside*nside)
    wind_i = np.zeros(12*nside*nside)
    
    for region in windows:
        d_ = ft.read(windows[i])
        wind_i[d_['hpix']] += d_['weight'] 
        count_i[d_['hpix']] += 1.0

    
    output_path = f'/global/cscratch1/sd/arosado/test_sysnet/prepared_contaminated_mocks/windows_clean{t}/nnwindow_{i}.hp{nside}.fits'
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
         os.makedirs(output_dir)
    
    is_good = count_i > 0.0
    wind_i[is_good] = wind_i[is_good] / count_i[is_good]
    wind_i[~is_good] = hp.UNSEEN
    print(f'wrote {output_path}')
    hp.write_map(output_path, wind_i, dtype=np.float64, fits_IDL=False, overwrite=True)    
