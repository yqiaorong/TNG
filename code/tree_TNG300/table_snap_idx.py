"""We had the mass table but we didn't know 
1. the corresponding halo index at that snapshot;
2. The linking subhalo masses.
This scripts aim to give 
1. the halo index at each snapshot which matches the mass table;
2, the subhalo masses."""

import os
import numpy as np
from tqdm import tqdm
import pandas as pd
import illustris_python as il

def get_feature_at_snapX(simpath, df, field):  
    """This function searches the target feature of the subhalos with global indidces.
    e.g. SubhaloMass, 
         SubhaloGrNr (the parent halo's local index in snap X)
    """
    for irow, row in tqdm(df.iterrows(), desc='snaps'):

        snapnum = int(irow[5:])
        print(irow)
        
        # Load subhalo global index 
        Subhalo = il.groupcat.loadSubhalos(simpath+'output', snapnum, fields=field)
        print(f'{field} loaded')
        # -------------------------------------------------------------------------------------
        # Substitute the subhalo global index with the parent halo index in the snapshot X
        df.loc[irow] = [Subhalo[int(idx)] if pd.notna(idx) else np.nan for idx in row]
 
        # -------------------------------------------------------------------------------------
        print(f'{field} assigned')
    return df

boxsize, res = 205, 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
sim_type = 'Hydro'

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
# Make a copy of tot_df to save the subhalo mass
tot_df_subhalo = tot_df.copy()

# Now change the entries of these subhalo global index to their "halo index in SnapX'
tot_df = get_feature_at_snapX(basePath, tot_df, 'SubhaloGrNr')

# Now change the entries of these subhalo global index to their "subhalo masses'
tot_df_subhalo = get_feature_at_snapX(basePath, tot_df_subhalo, 'SubhaloMass')

# Save the dataframe
save_mass_dir = f'result/DMhalo_table_mass/TNG300/sim_{boxsize}_{res}_{sim_type}/'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)
    
final_index_values = ['snap_' + str(i) for i in [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]] 
tot_df.to_csv(f'{save_mass_dir}/halo_snap_idx_table.csv', index=final_index_values)
tot_df_subhalo.to_csv(f'{save_mass_dir}/subhalo_mass_table.csv', index=final_index_values)