"""All the data computed in this script is saved in result/DMhalo_density_profiles"""

import os
import argparse
from tqdm import tqdm
import pandas as pd
import numpy as np
from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default='DM-Arepo/MTNG-L500-4320-A', type=str)
parser.add_argument('--snapnum',default=129, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add formation time to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

snap_list = [51, 69, 94, 129, 151, 179, 214, 237, 264]

# Load the index table
# ------------------------------------------------------------------------------------------------
FPGr_fnames = [f"snap_{snap}_FPGr.npy" for snap in snap_list]  
idx_table_paths = [f"result/DMhalo_mass_table_new/{args.sim}/{fname}" for fname in FPGr_fnames]  

# Load index tables into a dictionary
df_idx = {}
for snap, path in tqdm(list(zip(snap_list, idx_table_paths)), desc='load idx table'):
    df_idx [f"snap_{snap}"] = np.load(path) 


# Load the mass data
# ------------------------------------------------------------------------------------------------
FPGrMass_fnames = [f"snap_{snap}_FPGrMass.npy" for snap in snap_list]  
mass_table_paths = [f"result/DMhalo_mass_table_new/{args.sim}/{fname}" for fname in FPGrMass_fnames]  

# Load index tables into a dictionary
df_mass = {}
for snap, path in tqdm(list(zip(snap_list, mass_table_paths)), desc='load mass table'):
    df_mass[f"snap_{snap}"] = np.load(path)  


# Add the formation time to the snap data
# ------------------------------------------------------------------------------------------------
add_formation_time(args.sim, args.snapnum, df_idx, df_mass)