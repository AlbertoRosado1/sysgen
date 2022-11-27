import numpy as np

import os
import argparse
import sys
sys.path.append("/global/homes/a/arosado/sysgen/")
from utils import normalize_selection_func, make_hp, radec2hpix
from astropy.table import Table

parser = argparse.ArgumentParser() #AJRM
parser.add_argument("--zbin",help="zbin used: 0, 1, or 2.",default=0,type=int)
parser.add_argument("--phstart", help="mock realization to start from.",default=0,type=int)
parser.add_argument("--phend", help="mock realization to finish on.",default=24,type=int)
parser.add_argument("--mockdir", help="directory where contaminated mocks are",default=None)
args = parser.parse_args()

#import sys
zbin=args.zbin # int(sys.argv[1])
nside=256
print(zbin)
print('pycorr')

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
cont_dir = '/global/cscratch1/sd/arosado/test_sysnet/contaminated_mocks/' # AJRM

if args.mockdir == None:
    mock_dir = cont_dir
else:
    mock_dir = args.mockdir
print(f"using mocks from {mock_dir}")

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

    print('Randoms loaded: '+str(n+1)+' out of '+str(len(S_list)))
    print(len(ra_rand))

    n+=1

print('')
print('Total of randoms:')
print(len(ra_rand))

# Let's import pycorr

from pycorr import TwoPointCorrelationFunction, project_to_multipoles, utils, TwoPointCounter

# Let's import cosmoprimo

import cosmoprimo

#Defining cosmology

cosmo=cosmoprimo.fiducial.DESI()
cosmo.set_engine('camb')

#Converting coordinates

dists_rand=cosmo.comoving_radial_distance(z_rand)
randoms_positions=utils.sky_to_cartesian(np.array([ra_rand,dec_rand,dists_rand]))






for n in np.arange(args.phstart,args.phend+1): # AJRM added start en stop
    
    ph=str(n).zfill(3)
    
    # Let's load the data
    
    filename_orig='cutsky_LRG_z0.800_AbacusSummit_base_c000_ph'+ph
    filename=filename_orig+'.fits'
    #if not os.path.exists('results_pycorr/Xi_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy'):
    if not os.path.exists(mock_dir+'results_pycorr_corrected/Xi_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy'): #AJRM
    
        ra,dec,z,status= read(mock_dir+filename) # read(loc+filename)

        #a=(z>zmin) & (z<zmax) & (status & 2**1>0) & (status & 2**3>0) # Y5, main
        a=(z>zmin) & (z<zmax) & (status & 2**1>0) # AJRM raw Y5 slection 

        ra=ra[a]
        dec=dec[a]
        z=z[a]

        print('Data loaded')

        print('')
        print('Total of galaxies:')
        print(len(ra))
        
        ########################################################################################################
        # AJRM get predicted galaxy counts from NN, and convert to hpix map
        base_dir = os.path.join(os.getenv('CSCRATCH'), 'test_sysnet')
        cont_mocks_tables_dir = os.path.join(base_dir, 'contaminated_mocks_tables')
        mdir = os.path.join(cont_mocks_tables_dir, f'cutsky_LRG_z0.800_AbacusSummit_base_c000_ph{ph}')
        nn_table = Table.read(f'{mdir}/nn-weights.fits')
        npred = make_hp(np.mean(nn_table['weight'],axis=1), nn_table['hpix'], nside) # take mean along snapshots
        nn_normed = normalize_selection_func(npred,norm_method='mean')
        
        # AJRM select pix for nn-weights from ra and dec of contaminated mock
        hpix_for_nn = radec2hpix(nside, ra, dec)
        nn = nn_normed[hpix_for_nn]
        nn[nn==0] = 1 # AJRM If there are no predicted galaxy counts from NN then just leave alone 
        data_weights = 1/nn # AJRM contamination is multplicative so we divide by predicted galaxy counts
        ########################################################################################################
        
        #Converting coordinates

        dists=cosmo.comoving_radial_distance(z)
        data_positions=utils.sky_to_cartesian(np.array([ra,dec,dists]))

        edges=(np.linspace(0,200,201),np.linspace(-1,1,241)) # bin edges


        #fn=f'{mock_dir}results_pycorr/R1R2'
        fn=f'{mock_dir}results_pycorr_corrected/R1R2' # AJRM
        fn+='_zmin'+str(zmin)+'_zmax'+str(zmax)

        import os

        if os.path.exists(fn+'.npy'):

            print('R1R2 already exists, so no need to compute it again')

            R1R2=TwoPointCounter.load(fn+'.npy')

            result=TwoPointCorrelationFunction('smu',edges,data_positions1=data_positions,data_weights1=data_weights,
                                               randoms_positions1=randoms_positions,engine='corrfunc',nthreads=64,R1R2=R1R2,mpicomm=comm,dtype='f8')
        else:

            result=TwoPointCorrelationFunction('smu',edges,data_positions1=data_positions,data_weights1=data_weights,
                                               randoms_positions1=randoms_positions,engine='corrfunc',nthreads=64,mpicomm=comm,dtype='f8')

            result.R1R2.save(fn)

        # Let us project to multipoles (monopole, quadruple, hexadecapole)
        ells=(0,2,4)
        s,xiell=project_to_multipoles(result,ells=ells)
        
        #result.save('results_pycorr/Xi_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy')
        result.save(mock_dir+'results_pycorr_corrected/Xi_'+filename_orig+'_zmin'+str(zmin)+'_zmax'+str(zmax)+'.npy') #AJRM

        print('')
        print(n)