"""The answers are 13, 25, 40, 67, 99"""

import h5py
import illustris_python as il
import numpy as np
import argparse
from tqdm import tqdm
from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM', default=None, type=str)
args = parser.parse_args()

print('')
print(f'>>> Calculate cosmological time  <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



boxsize, res = 205, 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if args.DM == 'DM':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'_DM/output/'
elif args.DM == 'Hydro':
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'/output/'


# Load redshift
snap_list = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
redshifts = []
for snap in tqdm(snap_list, desc='load z'):
    with h5py.File(il.snapshot.snapPath(basePath+'output', snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z = 1 / scale_factor - 1
        redshifts.append(z)
redshifts = np.array(redshifts)

# Calculate t_dyn
prev_snap_idx, prev_z = calc_tdyn(redshifts, basePath, 13)
print(snap_list[prev_snap_idx], redshifts[prev_snap_idx], prev_z)