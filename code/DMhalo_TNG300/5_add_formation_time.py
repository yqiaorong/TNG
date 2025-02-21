"""All the data computed in this script is saved in result/DMhalo_density_profiles"""

import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from accret_func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default='TNG300/sim_205_1250_DM/', type=str)
parser.add_argument('--snapnum',default=25, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add formation time to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load the index table
# ------------------------------------------------------------------------------------------------
idx_table_dir = f'result/DMhalo_table_mass/{args.sim}/snap_idx_table.csv'
df_idx = pd.read_csv(idx_table_dir, index_col=0)
# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)

# Remove column names and convert df_idx to dictionary
idx_dict = {k: v.to_numpy() for k, v in tqdm(df_idx.iterrows())}
print('idx dict ready')

# Load the mass data
# ------------------------------------------------------------------------------------------------
mass_table_dir = f'result/DMhalo_table_mass/{args.sim}/mass_table.csv'
df_mass = pd.read_csv(mass_table_dir, index_col=0)

# Sort the mass so it matches the column of the index table
df_mass_sort = df_mass[df_idx.columns]
# Check if all of two columns are equal
print(np.all(df_mass_sort.columns == df_idx.columns))
del df_mass, df_idx 

# Remove column names and convert df_mass_sort to dictionary
mass_dict = {k: v.to_numpy() for k, v in tqdm(df_mass_sort.iterrows())}
del df_mass_sort
print('mass dict ready')

# Add the formation time to the snap data
# ------------------------------------------------------------------------------------------------
add_formation_time(args.sim, args.snapnum, idx_dict, mass_dict)