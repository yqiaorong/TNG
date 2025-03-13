"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from accret_func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type',default=None, type=str)
parser.add_argument('--snapnum', default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add formation time to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

table_dir = f'result/DMhalo_table_mass/TNG300/sim_205_1250_{args.sim_type}/'

# ------------------------------------------------------------------------------------------------
# Load the halo masses snap index table
df_idx = pd.read_csv(table_dir+'halo_snap_idx_table.csv', index_col=0)

# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)

# Remove column names and convert df_idx to dictionary
idx_dict = {k: v.to_numpy() for k, v in tqdm(df_idx.iterrows())}
print('idx dict ready')


# ------------------------------------------------------------------------------------------------
# Load the mass data
df_mass = pd.read_csv(table_dir+'halo_mass_table.csv', index_col=0)
if np.all(df_mass.columns == df_idx.columns) == False:
	exit()
else:
	# Remove column names and convert df_mass_sort to dictionary
	mass_dict = {k: v.to_numpy() for k, v in tqdm(df_mass.iterrows())}
	del df_mass
	print('mass dict ready')


# ------------------------------------------------------------------------------------------------
# Load the linking subhalo mass table
df_subhalo_mass = pd.read_csv(table_dir+'subhalo_mass_table.csv', index_col=0)
if np.all(df_subhalo_mass.columns == df_idx.columns) == False:
    exit()
else:
	# Remove column names and convert df_subhalo_mass to dictionary
	subhalo_mass_dict = {k: v.to_numpy() for k, v in tqdm(df_subhalo_mass.iterrows())}
	del df_subhalo_mass
	print('subhalo mass dict ready')


# Add the formation time to the snap data
# ------------------------------------------------------------------------------------------------
add_formation_time(args, idx_dict, mass_dict, subhalo_mass_dict)