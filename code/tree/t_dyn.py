"""The answers are ~51, ~69, 94, 151, 214, 264"""

import h5py
import illustris_python as il
import numpy as np
import argparse
from tqdm import tqdm
from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim', default='DM', type=str)
args = parser.parse_args()

print('')
print(f'>>> Calculate cosmological time <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

simpath = f'{args.sim}-Arepo/MTNG-L500-4320-A/'
basePath = f'/virgotng/mpa/MTNG/{simpath}/'

# Load redshift
snap_list = [51, 69, 80, 94, 129, 151, 179, 214, 237, 264]
redshifts = []
for snap in tqdm(snap_list, desc='load z'):
    with h5py.File(il.snapshot.snapPath(basePath+'output', snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        redshifts.append(z)
redshifts = np.array(redshifts)

# Calculate t_dyn
prev_snap_idx, prev_z = calc_tdyn(redshifts, basePath, 69)
print(snap_list[prev_snap_idx], redshifts[prev_snap_idx], prev_z)