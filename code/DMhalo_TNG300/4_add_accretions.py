"""All the data computed in this script is saved in result/DMhalo_density_profiles"""

import pandas as pd
import numpy as np
from accret_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type', default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add accretion rate to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# ------------------------------------------------------------------------------------------------
# Load the halo masses snap index table
idx_table_dir = f'result/DMhalo_table_mass/TNG300/sim_205_1250_{args.sim_type}/snap_idx_table.csv'
df_idx = pd.read_csv(idx_table_dir, index_col=0)

# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)

# ------------------------------------------------------------------------------------------------
# Load the accretion rate data
accret_table_dir = f'result/DMhalo_table_mass/TNG300/sim_205_1250_{args.sim_type}/accretion_table_new.csv'
df_accret = pd.read_csv(accret_table_dir, index_col=0)

# Sort the accretion rate so it matches the column of the index table
df_accret_sort = df_accret[df_idx.columns]
# Check if all of two columns are equal
print(np.all(df_accret_sort.columns == df_idx.columns))
del df_accret

# Get the row table data
snap_idx_table = df_idx.loc[f'snap_{args.snapnum}']
snap_accret_table = df_accret_sort.loc[f'snap_{args.snapnum}']

# Add the accretion rate to the snap data
_ = add_accret(args.sim_type, args.snapnum, snap_idx_table, snap_accret_table,)