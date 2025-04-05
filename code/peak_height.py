"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import os
import numpy as np
from DMhalo_TNG300.func import *
from tqdm import tqdm
import argparse
from colossus.cosmology import cosmology
from colossus.lss import peaks

cosmology.setCosmology('planck15')

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='TNG300/sim_205_1250_Hydro/', type=str) # [TNG300/MTNG]
parser.add_argument('--snapnum', default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add peak height to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    z = data['z']
    halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 MSun}
    
    # Caculate the peak height
    peakHeight = peaks.peakHeight(halo_M_Mean200*10**10, z)
    data['peakHeight'] = peakHeight
    print(data['halo_M_Mean200'].shape, data['peakHeight'].shape)
    
    # Save the data
    print(data.keys())
    np.save(os.path.join(halos_dir, fname), data)