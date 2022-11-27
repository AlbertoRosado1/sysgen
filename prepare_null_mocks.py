# prepare contaminated mocks in tables for NN mitigation
import sys 
sys.path.append("/global/homes/a/arosado/")
from prepare_data import prepare_LRGmock_data, maps_dr9, hpixsum
from astropy.table import Table, vstack

import os
import utils as ut
import numpy as np
import fitsio as ft
from glob import glob

norm_method = 'mean'
downsampling= 'frac'
nside = 256
#ph = 0
phs = [i for i in range(25)]
hp_order = 'ring'
zmin, zmax = 0., 6.
selection = dict(main=0, nz=0, Y5=1, sv3=0)
number_random_files = 5 # 20

scratch = os.getenv('CSCRATCH')
null_dir = os.path.join(scratch,'test_sysnet','null_mocks')
tables_dir = os.path.join(scratch, 'test_sysnet','contaminated_mocks_tables')

# read systematics
print('concatenating systematics')
systematics = vstack([Table.read(f'{scratch}/rands/rongpu_imaging_maps/pixmap_{r}_nside_256_minobs_1_maskbits_1111213.fits', format='fits') for r in ['north','south']])

# collect randoms
#print(f"concatenating {number_random_files} randoms")
#randoms_fn = glob('/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/LRG/z0.800/cutsky_LRG_random_*.fits')
#randoms_concat = vstack([Table.read(rand_fn, format='fits') for rand_fn in randoms_fn[:number_random_files]])
fn = f"/global/cscratch1/sd/arosado/test_sysnet/cutsky_LRG_random_{number_random_files}_concat.fits"
randoms = ut.apply_mock_mask(fn, **selection)

for ph in phs:
    prep_data_dir = os.path.join(tables_dir,f"cutsky_LRG_z0.800_AbacusSummit_base_c000_ph{ph:003d}")
    
    if not os.path.exists(prep_data_dir):
        os.mkdir(prep_data_dir)
        print('made '+prep_data_dir)
        
    # get contaminated mock (contaminated mock by downsampling)
    mock_fn = os.path.join(null_dir, f'cutsky_LRG_z0.800_AbacusSummit_base_c000_ph{ph:003d}.fits')
    mock = Table.read(mock_fn)
    print(f"preparing {mock_fn}")

    # prepare data for NN
    prep_data = prepare_LRGmock_data(mock, randoms, systematics, hp_order=hp_order, zmin=zmin, zmax=zmax, nside=nside, columns=maps_dr9)
    #fn = os.path.join(scratch, 'test_sysnet', 'prepared_contaminated_mocks', f'cutsky_LRG_z0.800_AbacusSummit_base_c000_ph{ph:003d}.fits')
    fn = os.path.join(prep_data_dir, f'null_sh_table_{nside}.fits')
    ft.write(fn, prep_data)
    print(f'prepped data and saved {fn}')