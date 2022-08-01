import os
import healpy as hp
import numpy as np
import fitsio as ft
from glob import glob
from time import time
from astropy.table import Table, vstack
import sys
sys.path.append('/global/homes/a/arosado/LSSutils/')
from lssutils.stats.cl import get_cl
import lssutils.utils as ut

# DEFINE FUNCTIONS
def radec2hpix(nside, ra, dec):
    """ change ra,dec to HEALPix hpix with ring ordering """
    theta, phi = np.deg2rad(90.-dec), np.deg2rad(ra)
    hpix = hp.ang2pix(nside, theta, phi, nest=False)     
    return hpix

def project2hp(nside, mock, weight=None, return_hpix=False):
    """ count # of objects in HEALPix """
    ra = mock['RA']
    dec = mock['DEC']
    hpix = radec2hpix(nside, ra, dec)
    hpmap = np.bincount(hpix, weights=weight, minlength=12*nside*nside)
    if return_hpix:
        return hpmap, hpix
    else:
        return hpmap
    
def read_data(data, nside=256):
    d_ = data
    ng = ut.make_hp(nside, d_['hpix'], d_['label'])
    nr = ut.make_hp(nside, d_['hpix'], d_['fracgood'])
    mask = ut.make_hp(nside, d_['hpix'], 1.0) > 0.5
    syst = ut.make_hp(nside, d_['hpix'], np.log(d_['features'][:, 0]))
    return ng, nr, mask, syst
    
def normalize_selection_func(ngal_pred):
    # get nside from the selection function
    nside = hp.get_nside(ngal_pred)
    
    # normalize the selection function to [0, 1]
    good = ngal_pred>0
    vmin, vmax = np.percentile(ngal_pred[good], [0, 100])
    selection_func = np.zeros_like(ngal_pred)
    selection_func[good] = (ngal_pred[good]-vmin) / (vmax-vmin)
    return selection_func

def downsample(selection_func, mock):
    """ downsample a mock catalog with a given selection function """
    nside = hp.get_nside(selection_func)
    hpix = radec2hpix(nside, mock['RA'], mock['DEC'])
    prob = selection_func[hpix]    
    good = np.random.uniform(size=mock.size) < prob
    return mock[good]

def _desi_mock_filename(tracer='LRG', ph=0):
    """Collect the name of DESI Mocks in NERSC."""
    mock_dir='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'
    z_dic={'LRG':0.8,'ELG':1.1,'QSO':1.4}
    fname=f'{mock_dir}{tracer}/z{z_dic[tracer]:5.3f}/cutsky_{tracer}_z{z_dic[tracer]:5.3f}_AbacusSummit_base_c000_ph{ph:003d}.fits'
    return fname

def _selection_function_filename(tracer='LRG', ind=0):
    """Collect selection function obtained from NN."""
    tracer = tracer.lower() # after testing make a directory by tracer
    fname=f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean/nnwindow_{ind}.hp256.fits'
    return fname

def make_hp(hpix, value, nside):
    """ A Function to create a HEALPix map
    """
    m_ = np.zeros(12*nside*nside)
    m_[:] = np.nan
    m_[hpix] = value
    
    return m_


# BEGIN CODE
output_dir = f'/global/cscratch1/sd/arosado/test_sysnet/results/cell'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# collect real data used to train NN
data = glob(f'/global/cscratch1/sd/arosado/nlrg_features_*_256.fits')
# stack all regions:
tables = []
for d_fn in data:
    d = Table.read(d_fn)
    #print(len(d))
    tables.append(d)
tables = vstack(tables) 

# read the mock catalog
mock_name = _desi_mock_filename(tracer='LRG', ph=0)
mock = ft.read(mock_name)

# collect randoms
randoms_path = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/LRG/z0.800/cutsky_LRG_random_S1000_1X.fits'
randoms = ft.read(randoms_path)
# keep only randoms with correct n(z) distribution (need for power spectrum measurement)
randoms = randoms[randoms['STATUS'] & 2**0 != 0]

nwindows = 1000 # total number of windows to pick from
only_choose = 100 # number of windows to use

np.random.seed(0)
pick_selfn_randomly = np.random.permutation(np.arange(nwindows)) # randomize selection functions to use for contamination

t0 = time()
for i, ind in enumerate(pick_selfn_randomly[:only_choose]):
    print(f'running {i+1}/{only_choose}. time elapsed: {time() - t0} s')
    # read selection function (obtained from real data)
    sel_fn = _selection_function_filename(tracer='LRG', ind=ind)
    ngal_pred = hp.read_map(sel_fn)
    
    # get nside from the selection function
    nside = hp.get_nside(ngal_pred)
    
    # normalize the selection function to [0, 1]
    selection_func = normalize_selection_func(ngal_pred)

    # subsample the mock catalog, and project to HEALPix
    mock_after = downsample(selection_func, mock)
    mock_after_hmap, hpix = project2hp(nside, mock_after, return_hpix=True)
    mask = make_hp(hpix , 1.0, nside)
    mask = mask > 0.5  # to make it binary
    
    # project randoms to healpix map
    randoms_hmap = project2hp(nside, randoms)
    cl_ = get_cl(mock_after_hmap, randoms_hmap, mask)
    
    # save results
    output_fn = f'{output_dir}/cell_{ind}.hp256.npy'
    print(f'saving {output_fn}')
    np.save(output_fn, cl_)
    
    
print(f'finished at {time() - t0} s')   
