import os
import numpy as np
import healpy as hp
from glob import glob
import fitsio as ft

# Run this script ONLY after running combine_nn_windows.py .
# nnwindow_*.hp{nside}.fits should be readable by hp.read_map, i.e. they should be healpix maps

t = '_test1'
nside = 256
windows = glob(f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean{t}/nnwindow_*.hp{nside}.fits')

count = np.zeros(12*nside*nside)
wind = np.zeros(12*nside*nside) 

for window in windows:
    w = hp.read_map(window)
    is_seen = (w != hp.UNSEEN)
    wind[is_seen] += w[is_seen]  
    count[is_seen] += 1.0

output_path = f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean{t}/nn-weights.hp{nside}.fits'
output_dir = os.path.dirname(output_path)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
is_good = count > 0.0
wind[is_good] = wind[is_good] / count[is_good]
wind[~is_good] = hp.UNSEEN
print(f'wrote {output_path}')
hp.write_map(output_path, wind, dtype=np.float64, fits_IDL=False, overwrite=True)    
