"""Add data to already compiled DMhalo data."""

import os
import argparse
import numpy as np
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default='Hydro-Arepo/MTNG-L500-4320-A',type=str)
parser.add_argument('--snapnum',      default=129,                      type=int)
parser.add_argument('--bin_start',    default=3.5,                     type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=4,                     type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

print('')
print('>>> Add data to compile DM halo density profiles <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



### Load h, scale_factor and z ###
load_dir = f'/nfs/mvogelsblab001/Users/s_qyu/{args.save_root_dir}/{args.sim}/snap_{args.snapnum}/final_densities/'
data = np.load(f'{load_dir}/bin-30-35.npy', allow_pickle=True).item()


    
### Add data to final compiled data
compiled_data = np.load(f'{load_dir}/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}.npy',
                        allow_pickle=True).item()

compiled_data['h'] = data['h']
compiled_data['scale_factor'] = data['scale_factor']
compiled_data['z'] = data['z']

np.save(f'{load_dir}/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}', compiled_data)
print('saved')