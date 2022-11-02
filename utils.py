import os
import fitsio as ft
from astropy.table import Table
import numpy as np
import healpy as hp
from scipy import interpolate


def radec2hpix(nside, ra, dec):
    """ change ra,dec to HEALPix hpix with ring ordering """
    theta, phi = np.deg2rad(90.-dec), np.deg2rad(ra)
    hpix = hp.ang2pix(nside, theta, phi, nest=False)     
    return hpix

def hpix2radec(nside, hpix):
    """
    Function transforms HEALPix index (ring) to RA, DEC
    
    parameters 
    ----------
    nside : int
        
    hpix : array_like
        HEALPix indices
    
    returns
    -------
    ra : array_like
        right ascention in deg
        
    dec : array_like
        declination in deg
        
    """
    theta, phi = hp.pixelfunc.pix2ang(nside, hpix)
    return np.degrees(phi), 90-np.degrees(theta)

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
    
def normalize_selection_func(ngal_pred, norm_method='mean'):
    assert norm_method in 'mean_percentile', "normalization method must be 'mean' or 'percentile'"
    # get nside from the selection function
    nside = hp.get_nside(ngal_pred)
    
    good = ngal_pred>0
    # normalize the selection function to [0, 1]
    if norm_method=='mean':
        print('contamination using mean')
        vmin, vmax = ngal_pred[good].min(), ngal_pred[good].mean() #np.percentile(ngal_pred[good], [0, 100])
        selection_func = np.zeros_like(ngal_pred)
        selection_func[good] = (ngal_pred[good]/vmax)#-vmin) / (vmax-vmin)
    else:
        print('contamination using percentile')
        vmin, vmax = np.percentile(ngal_pred[good], [0, 100])
        selection_func = np.zeros_like(ngal_pred)
        selection_func[good] = (ngal_pred[good]-vmin) / (vmax-vmin)
    return selection_func

def downsample(selection_func, mock, downsampling='mean', seed=None):
    """ downsample a mock catalog with a given selection function """
    assert downsampling in 'mean_frac', "downsampling method must be 'mean' or 'frac'"
    nside = hp.get_nside(selection_func)
    hpix = radec2hpix(nside, mock['RA'], mock['DEC'])
    prob = selection_func[hpix]
    if downsampling=='frac':
        print('using frac when downsampling')
        raw = Table.read('LRG_nz_raw.fits')
        main = Table.read('LRG_nz_main.fits')
        fmain = interpolate.interp1d(main['zmid'], main['n(z)'], fill_value='extrapolate')
        fraw = interpolate.interp1d(raw['zmid'], raw['n(z)'], fill_value='extrapolate')
        dndz_main = fmain(mock['Z'])
        dndz_raw = fraw(mock['Z'])
        frac = dndz_main/dndz_raw
        prob *= frac
    rng = np.random.seed(seed=seed)
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

def _selection_function_filename_mean(tracer='LRG'):
    """Collect selection function obtained from NN."""
    tracer = tracer.lower() # after testing make a directory by tracer
    fname=f'/global/cscratch1/sd/arosado/test_sysnet/windows_clean_test1/nn-weights.hp256.fits'
    return fname

def make_hp(value, hpix, nside, fill_with=np.nan):
    """ A Function to create a HEALPix map
    """
    m_ = np.zeros(12*nside*nside)
    m_[:] = fill_with
    m_[hpix] = value
    
    return m_

def get_mock_hpmap(contaminated=False, selection_fn=None, nside=256, tracer='LRG', ph=0, main=0, nz=1, Y5=1, sv3=0):
    
    # read the mock catalog
    mock_name = _desi_mock_filename(tracer=tracer, ph=ph)
    print(f'using {mock_name}')
    mock = apply_mock_mask(mock_name, main=main, nz=nz, Y5=Y5, sv3=sv3) # mask Y5 by default
    
    if contaminated:
        assert selection_fn is not None, "provide selection function"
        # normalize the selection function to [0, 1]
        selection_func = normalize_selection_func(selection_fn)
        # subsample the mock catalog, and project to HEALPix
        mock = downsample(selection_func, mock)
        
    mock_hpmap = project2hp(nside, mock)
    return mock_hpmap

def get_mock(contaminated=False, selection_fn=None, norm_method='mean', downsampling='mean', nside=256, seed=42,
             tracer='LRG', ph=0, return_hpix=False, main=0, nz=0, Y5=1, sv3=0):
    if selection_fn is not None:
        nside = hp.get_nside(selection_fn)
    # read the mock catalog
    mock_name = _desi_mock_filename(tracer=tracer, ph=ph)
    print(f'using {mock_name}')
    mock = apply_mock_mask(mock_name, main=main, nz=nz, Y5=Y5, sv3=sv3) # mask Y5 by default
    
    if contaminated:
        assert selection_fn is not None, "provide selection function"
        print('contaminating mock')
        # normalize the selection function to [0, 1]
        selection_func = normalize_selection_func(selection_fn, norm_method=norm_method)
        # subsample the mock catalog, and project to HEALPix
        mock = downsample(selection_func, mock,  downsampling=downsampling, seed=seed)
        
    if return_hpix:
        hpix = radec2hpix(nside, mock['RA'], mock['DEC'])
        return mock, hpix
    else:
        return mock
    
def mock_mask(main=0, nz=0, Y5=0, sv3=0):
    #https://desi.lbl.gov/trac/attachment/wiki/CosmoSimsWG/FirstGenerationMocks/STATUS_script_example.py
    return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)

def apply_mock_mask(filename, main=1, nz=0, Y5=1, sv3=0):
    #https://desi.lbl.gov/trac/attachment/wiki/CosmoSimsWG/FirstGenerationMocks/STATUS_script_example.py
    ## Default value chooses the Y5 footprint and downsampled to input n(z)
    data = ft.read(filename) #Table.read(filename)
    STATUS = data["STATUS"]

    ### Array of indices
    idx = np.arange(len(STATUS))
    
    #apply mask to get desired sample
    print(f'applying mask: main={main}, nz={nz}, Y5={Y5}, sv3={sv3}')
    mask = mock_mask(main=main, nz=nz, Y5=Y5, sv3=sv3) # mask Y5 by default
    idx_good = idx[(STATUS & (mask))==mask]
    return data[idx_good]