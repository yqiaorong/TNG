import pandas as pd
import numpy as np
from accret_func import *

# Initial params
sim_type = 'DM'
snapnum = 99

# ------------------------------------------------------------------------------------------------
# Load the halo masses snap index table
idx_table_dir = f'result/DMhalo_mass_table/TNG300/sim_205_1250_{sim_type}/snap_idx_table.csv'
df_idx = pd.read_csv(idx_table_dir, index_col=0)

# Convert non nan entries in df_idx to int
df_idx = df_idx.fillna(-1)
df_idx = df_idx.astype(int) 
df_idx = df_idx.replace(-1, np.nan)

# ------------------------------------------------------------------------------------------------
# # Load the halo masses in the mass table at the same redshift
# mass_table_dir = f'result/DMhalo_mass_table/TNG300/sim_205_1250_{sim_type}/mass_table.csv'
# df_mass = pd.read_csv(mass_table_dir, index_col=0)

# # Sort the accretion rate so it matches the column of the index table
# df_mass_sort = df_mass[df_idx.columns]
# # Check if all of two columns are equal
# print(np.all(df_mass_sort.columns == df_idx.columns))
# del df_mass

# ------------------------------------------------------------------------------------------------
# Load the accretion rate data
accret_table_dir = f'result/DMhalo_mass_table/TNG300/sim_205_1250_{sim_type}/accretion_table.csv'
df_accret = pd.read_csv(accret_table_dir, index_col=0)

# Sort the accretion rate so it matches the column of the index table
df_accret_sort = df_accret[df_idx.columns]
# Check if all of two columns are equal
print(np.all(df_accret_sort.columns == df_idx.columns))
del df_accret

# Get the row table data
snap_idx_table = df_idx.loc[f'snap_{snapnum}']
snap_accret_table = df_accret_sort.loc[f'snap_{snapnum}']

# Add the accretion rate to the snap data
_ = add_accret(sim_type, snapnum, snap_idx_table, snap_accret_table,)