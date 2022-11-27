import numpy as np

import os
import sys
import argparse
parser = argparse.ArgumentParser() #AJRM
parser.add_argument("--zbin",help="zbin used: 0, 1, or 2.",default=0,type=int)
parser.add_argument("--phstart", help="mock realization to start from.",default=0,type=int)
parser.add_argument("--phend", help="mock realization to finish on.",default=24,type=int)
parser.add_argument("--mockdir", help="directory where mocks are",default=None)
args = parser.parse_args()

#import sys
zbin= args.zbin # int(sys.argv[1])

# print(zbin)
# print('pypower')

from mpi4py import MPI

use_mpi=True

if use_mpi:
    comm=MPI.COMM_WORLD
else:
    comm=MPI.COMM_SELF

import fitsio

# Read input catalogs, scattering on all MPI ranks

def read(fn,columns=('RA','DEC','Z','STATUS'),ext=1,mpicomm=comm):
    gsize=fitsio.FITS(fn)[ext].get_nrows()
    start,stop=mpicomm.rank*gsize // mpicomm.size,(mpicomm.rank+1)*gsize // mpicomm.size
    tmp=fitsio.read(fn,ext=ext,columns=columns,rows=range(start,stop))
    return [tmp[col] for col in columns]

# Let's define zmin and zmax

if zbin==0:
    zmin=0.4
    zmax=0.6
elif zbin==1:
    zmin=0.6
    zmax=0.8
elif zbin==2:
    zmin=0.8
    zmax=1.1

# Location of the mocks 
loc='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/LRG/z0.800/' # location of original mocks, and randoms
null_dir = '/global/cscratch1/sd/arosado/test_sysnet/null_mocks/' # AJRM

if args.mockdir == None:
    mock_dir = null_dir
else:
    mock_dir = args.mockdir


# Let's load the randoms

NX_random=20

S_list=range(100,(1+NX_random)*100,100)

n=0

for S in S_list:
    
    S=str(S)
    
    filename_rand='cutsky_LRG_random_S'+S+'_1X.fits'
    
    if n==0:
        ra_rand,dec_rand,z_rand,status_rand=read(loc+filename_rand)

        #a_rand=(z_rand>zmin) & (z_rand<zmax) & (status_rand & 2**1>0) & (status_rand & 2**3>0)
        a_rand=(z_rand>zmin) & (z_rand<zmax) & (status_rand & 2**1>0) # AJRM raw Y5 slection

        ra_rand=ra_rand[a_rand]
        dec_rand=dec_rand[a_rand]
        z_rand=z_rand[a_rand]

    else:
        RA_RAND,DEC_RAND,Z_RAND,STATUS_RAND=read(loc+filename_rand)

        #A_RAND=(Z_RAND>zmin) & (Z_RAND<zmax) & (STATUS_RAND & 2**1>0) & (STATUS_RAND & 2**3>0)
        A_RAND=(Z_RAND>zmin) & (Z_RAND<zmax) & (STATUS_RAND & 2**1>0) # AJRM raw Y5 slection

        RA_RAND=RA_RAND[A_RAND]
        DEC_RAND=DEC_RAND[A_RAND]
        Z_RAND=Z_RAND[A_RAND]

        ra_rand=np.concatenate([ra_rand,RA_RAND])
        dec_rand=np.concatenate([dec_rand,DEC_RAND])
        z_rand=np.concatenate([z_rand,Z_RAND])

    n+=1

# Let's import pypower

from pypower import CatalogFFTPower, utils, setup_logging

# To activate logging

setup_logging()

# Let's import cosmoprimo

import cosmoprimo

#Defining cosmology

cosmo=cosmoprimo.fiducial.DESI()
cosmo.set_engine('camb')

#Converting coordinates

dists_rand=cosmo.comoving_radial_distance(z_rand)
randoms_positions=utils.sky_to_cartesian(np.array([ra_rand,dec_rand,dists_rand]))




nmesh=1024

if zbin==2:
    boxsize=5000
else:
    boxsize=4000

edges={'min':0,'step':0.001}
ells=(0,2,4)

for n in np.arange(args.phstart,args.phend+1): # AJRM added start en stop
    
    ph=str(n).zfill(3)
    
    # Let's load the data
    
    filename_orig='cutsky_LRG_z0.800_AbacusSummit_base_c000_ph'+ph
    filename=filename_orig+'.fits'
    
    #if not os.path.exists('results_pypower_'+str(nmesh)+'/Pk_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy'):
    if not os.path.exists(mock_dir+'results_pypower_'+str(nmesh)+'/Pk_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy'):
        print(f"using mocks from {mock_dir+filename}")
        ra,dec,z,status=read(mock_dir+filename)

        #a=(z>zmin) & (z<zmax) & (status & 2**1>0) & (status & 2**3>0)
        a=(z>zmin) & (z<zmax) & (status & 2**1>0) # AJRM raw Y5 slection 

        ra=ra[a]
        dec=dec[a]
        z=z[a]

        #Converting coordinates

        dists=cosmo.comoving_radial_distance(z)
        data_positions=utils.sky_to_cartesian(np.array([ra,dec,dists]))

        result=CatalogFFTPower(data_positions1=data_positions,randoms_positions1=randoms_positions,
                               edges=edges,ells=ells,interlacing=2,boxsize=boxsize,nmesh=nmesh,resampler='tsc',
                               los=None,position_type='xyz',mpicomm=comm)

        #result.save('results_pypower_'+str(nmesh)+'/Pk_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy')
        result.save(mock_dir+'results_pypower_'+str(nmesh)+'/Pk_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy')