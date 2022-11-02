# concatenate randoms for mocks
from astropy.table import Table, vstack
import os
import utils as ut
import fitsio as ft
from glob import glob

number_random_files = 5

# collect randoms
print(f"concatenating {number_random_files} randoms")
randoms_fn = glob('/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/LRG/z0.800/cutsky_LRG_random_*.fits')

randoms_concat = Table()
for rand_fn in randoms_fn[:number_random_files]:
    print(f"concatenating {rand_fn}")
    randoms_concat = vstack([randoms_concat, Table.read(rand_fn)])
    
#randoms_concat = vstack([Table.read(rand_fn, format='fits') for rand_fn in randoms_fn[:number_random_files]])
out_fn = f"/global/cscratch1/sd/arosado/test_sysnet/cutsky_LRG_random_{number_random_files}_concat.fits"
randoms_concat.write(out_fn, format='fits')
print(f"created randoms file {out_fn}")