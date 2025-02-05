"""We had the mass table but we didn't know the corresponding halo index at that snapshot. 
This scripts aim to give the halo index at each snapshot which matches the mass table."""

import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import illustris_python as il
# from numba import njit, prange

# @njit(parallel=True)
# def assign_idx(row, SubhaloGrNr):
#     result = np.empty(len(row), dtype=np.float64) 
#     for i in prange(len(row)):  
#         if not np.isnan(row[i]):  
#             result[i] = SubhaloGrNr[int(row[i])]
#         else:
#             result[i] = np.nan 
#     return result

def get_halo_mass_idx_at_snapX(simpath, df):  
    
    for irow, row in tqdm(df.iterrows(), desc='snaps'):

        snapnum = int(irow[5:])
        print(irow)
        
        # Load subhalo global index 
        SubhaloGrNr = il.groupcat.loadSubhalos(simpath+'output', snapnum, fields='SubhaloGrNr')
        print('Index loaded')
        # -------------------------------------------------------------------------------------
        # Substitute the subhalo global index with the parent halo index in the snapshot X
        df.loc[irow] = [SubhaloGrNr[int(idx)] if pd.notna(idx) else np.nan for idx in row]
        # df.loc[irow] = assign_idx(row.to_numpy(), SubhaloGrNr)
        # -------------------------------------------------------------------------------------
        print('Index assigned')
    return df

boxsize, res = 205, 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
sim_type = 'DM'
if sim_type == 'DM':
    # DM
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'_DM/'
    idx_dir = f'result/DMhalo_table_chunk_idx/TNG300/sim_{boxsize}_{res}_DM/'
else:
    # hydro
    basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+'/'
    idx_dir = f'result/DMhalo_table_chunk_idx/TNG300/sim_{boxsize}_{res}_Hydro/'

# Concatenate all index df
df_list = os.listdir(idx_dir)
tot_df_list = []
for file in tqdm(df_list, desc='concatenate chunk df'):
    df = pd.read_csv(f'{idx_dir}/{file}', index_col=0)  # Assuming the first column is the index
    tot_df_list.append(df)
tot_df = pd.concat(tot_df_list, axis=1)

# Now change the entries of these global index to their masses
tot_df = get_halo_mass_idx_at_snapX(basePath, tot_df)

# Save the dataframe
save_mass_dir = f'result/DMhalo_table_mass/TNG300/sim_{boxsize}_{res}_{sim_type}/'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)
final_index_values = ['snap_' + str(i) for i in [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]] 
tot_df.to_csv(f'{save_mass_dir}/snap_idx_table.csv', index=final_index_values)