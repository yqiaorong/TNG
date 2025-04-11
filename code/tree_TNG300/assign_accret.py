"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import pandas as pd
import numpy as np
from accret_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type', default='Hydro', type=str)
parser.add_argument('--snapnum',  default=None,    type=int)
args = parser.parse_args()

print('')
print(f'>>> Add accretion rate to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


table_dir = f'result/DMhalo_table/TNG300/sim_205_1250_{args.sim_type}/'

# ------------------------------------------------------------------------------------------------
# Load the halo masses snap index table
df_idx = pd.read_csv(table_dir+'halo_index_table.csv', index_col=0)

# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)
print('Halo local index succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Load the accretion rate data
df_accret = pd.read_csv(table_dir+'accretion_table.csv', index_col=0) 
if np.all(df_accret.columns == df_idx.columns) == False:
    exit()
else:
    print('Accretion rate table succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Load the mass data
df_mass = pd.read_csv(table_dir+'halo_mass_table.csv', index_col=0) # Comoving mass
if np.all(df_mass.columns == df_idx.columns) == False:
	exit()
else:
    print('Halo mass table succesfully loaded')
    
# ------------------------------------------------------------------------------------------------
# Load the subhalo mass data
df_submass = pd.read_csv(table_dir+'subhalo_mass_table.csv', index_col=0) # Comoving mass
if np.all(df_submass.columns == df_idx.columns) == False:
	exit()
else:
    print('SUbhalo mass table succesfully loaded')

# ------------------------------------------------------------------------------------------------
# Get the row table data
snap_idx_table    = df_idx.loc[f'snap_{args.snapnum}']
snap_accret_table = df_accret.loc[f'snap_{args.snapnum}']
snap_mass_table   = df_mass.loc[f'snap_{args.snapnum}'] # Comoving mass!!!
snap_submass_table = df_submass.loc[f'snap_{args.snapnum}'] # Comoving mass!!!
del df_idx, df_accret, df_mass, df_submass

# ------------------------------------------------------------------------------------------------
# Add the accretion rate to the snap data
_ = add_accret(args.sim_type, args.snapnum, 
               snap_idx_table, snap_mass_table, snap_accret_table, snap_submass_table)