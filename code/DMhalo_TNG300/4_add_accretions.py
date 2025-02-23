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


table_dir = f'result/DMhalo_table_mass/TNG300/sim_205_1250_{args.sim_type}/'

# ------------------------------------------------------------------------------------------------
# Load the halo masses snap index table
df_idx = pd.read_csv(table_dir+'halo_snap_idx_table.csv', index_col=0)

# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)
print('Halo local index succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Load the linking subhalo mass table
df_subhalo_mass = pd.read_csv(table_dir+'subhalo_mass_table.csv', index_col=0)
if np.all(df_subhalo_mass.columns == df_idx.columns) == False:
    exit()
else:
    print('Subhalo mass table succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Load the accretion rate data
df_accret = pd.read_csv(table_dir+'accretion_table_new.csv', index_col=0)
if np.all(df_accret.columns == df_idx.columns) == False:
    exit()
else:
    print('Accretion rate table succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Load the mass data
df_mass = pd.read_csv(table_dir+'halo_mass_table.csv', index_col=0)
if np.all(df_mass.columns == df_idx.columns) == False:
	exit()
else:
    print('Mass table succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Get the row table data
snap_idx_table = df_idx.loc[f'snap_{args.snapnum}']
snap_accret_table = df_accret.loc[f'snap_{args.snapnum}']
snap_mass_table = df_mass.loc[f'snap_{args.snapnum}']
snap_subhalo_mass_table = df_subhalo_mass.loc[f'snap_{args.snapnum}']
del df_idx, df_accret, df_subhalo_mass, df_mass

# ------------------------------------------------------------------------------------------------
# Add the accretion rate to the snap data
_ = add_accret(args.sim_type, args.snapnum, 
               snap_idx_table, snap_mass_table, snap_accret_table, snap_subhalo_mass_table)