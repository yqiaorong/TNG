"""This script stores the subhalo snap index in a dataframe in result/DMhalo_table_chunk_idx/.
Then, it substitutes the subhalo snap index with the halo snap index and saves the dataframe in result/DMhalo_table/.
Finally, it substitutes the halo snap index with the halo mass and saves the dataframe in result/DMhalo_table/.
"""

import os
import argparse
from tqdm import tqdm
from func import *
import numpy as np
import pandas as pd
import illustris_python as il

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,   type=int)
parser.add_argument('--res',    default=1250,  type=int)
parser.add_argument('--sim_type', default='Hydro',type=str)
parser.add_argument('--bin_start',default=1,   type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',  default=4.5, type=float) # [10^{10+x} Msun/h]
args = parser.parse_args()

print('')
print(f'>>> Mass table <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if args.sim_type == 'DM':
    basePath = data_path + 'L%dn%dTNG'%(args.boxsize,args.res)+'_DM'
elif args.sim_type == 'Hydro':
    basePath = data_path + 'L%dn%dTNG'%(args.boxsize,args.res)
snapnum = 99 # The last snapshot, should not change
boxsize = args.boxsize
res     = args.res


# Select a subset of DM halos from groupcat
group_fields = ['Group_M_Mean200', 'GroupFirstSub']
Halos = il.groupcat.loadHalos(basePath+'/output', snapnum, fields=group_fields)

Group_M_Mean200 = Halos['Group_M_Mean200'] # [10^10 MSun/h]
GroupFirstSub = Halos['GroupFirstSub']
del Halos

# Using comoving mass to select subset
halo_global_idx = np.where((Group_M_Mean200 >= 10**args.bin_start) & 
                           (Group_M_Mean200 < 10**args.bin_end))[0]
num_halos = halo_global_idx.shape[0]
print(f'The number of halos in the given mass range: {num_halos}')


index_values = ['snap_' + str(i) for i in range(2, 99)]
final_index_values = ['snap_' + str(i) for i in [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]]  


# Save directory
save_chunk_dir = f'result/DMhalo_table_chunk_idx/TNG300/sim_{boxsize}_{res}_{args.sim_type}/'
if not os.path.exists(save_chunk_dir):
    os.makedirs(save_chunk_dir)

save_table_dir = f'result/DMhalo_table/TNG300/sim_{boxsize}_{res}_{args.sim_type}/'
if not os.path.exists(save_table_dir):
    os.makedirs(save_table_dir)
    
    
    
chunk_size = 100
chunk_indices = range(0, num_halos, chunk_size)



# from joblib import Parallel, delayed
# def process_chunk(chunk_idx):
#     start = chunk_idx
#     end = min(num_halos, chunk_idx + chunk_size)
#     df = pd.DataFrame(index=index_values)
    
#     for idx in halo_global_idx[start:end]:
#         subhalo_global_idx = GroupFirstSub[idx]
#         df = lifeline(basePath, subhalo_global_idx, df)
    
#     df = df.loc[final_index_values]
#     df.to_csv(f'{save_chunk_dir}/chunk_{chunk_idx}.csv', index=True)
#     return chunk_idx

# results = Parallel(n_jobs=-1)(delayed(process_chunk)(i) for i in tqdm(chunk_indices))



for chunk_idx in chunk_indices:
        
    ### Creat halo (first subhalo) lifeline dataframe ###
    df = pd.DataFrame(index=index_values)

    # Load treeX 
    for idx in tqdm(halo_global_idx[chunk_idx: min(num_halos, chunk_idx+chunk_size)], 
                    desc=f'make index df: chunk {chunk_idx}+'):
        subhalo_global_idx = GroupFirstSub[idx]
    
        # Get the dataframe of subhalo global index
        df = lifeline(basePath, subhalo_global_idx, df)

    # Only keep necessary snapshots
    df = df.loc[final_index_values]
    print(df.shape)

    # Save the dataframe
    df.to_csv(f'{save_chunk_dir}/chunk_{chunk_idx}.csv', index=True)



# Concatenate all dataframes
df_list = os.listdir(save_chunk_dir)
tot_df_list = []
for file in tqdm(df_list, desc='concatenate chunk df'):
    df = pd.read_csv(f'{save_chunk_dir}/{file}', index_col=0)  # Assuming the first column is the index
    tot_df_list.append(df)
tot_df = pd.concat(tot_df_list, axis=1)
print(tot_df)
print(tot_df.shape)
# ------------------------------------------------------------------------------------------------------
# tot_df is a dataframe of subhalo global index within snapshots.
# ------------------------------------------------------------------------------------------------------


# Now change the entries of these subhalo global index to their parent halo index
halo_index_df = get_halo_index_of_lifeline(basePath, tot_df.copy()) 
halo_index_df.to_csv(f'{save_table_dir}/halo_index_table.csv', index=final_index_values)
print(halo_index_df)
# -----------------------------------------------------------------------------------------------------
# halo_index_df is a dataframe of halo index within snapshots.
# -----------------------------------------------------------------------------------------------------


# Now change the entries of these halo index to their mass
halo_mass_df = get_group_field_of_lifeline(basePath, halo_index_df.copy(), group_field='Group_M_Mean200') # Comoving Mass!
halo_mass_df.to_csv(f'{save_table_dir}/halo_mass_table.csv', index=final_index_values)
print(halo_mass_df)
# -----------------------------------------------------------------------------------------------------
# halo_mass_df is a dataframe of halo mass within snapshots.
# -----------------------------------------------------------------------------------------------------


# Convert the subhalo global index to subhalo mass
subhalo_df = get_subhalo_field_of_lifeline(basePath, tot_df.copy(), subhalo_field='SubhaloMass')
subhalo_df.to_csv(f'{save_table_dir}/subhalo_mass_table.csv', index=final_index_values)
print(subhalo_df)
# -----------------------------------------------------------------------------------------------------
# subhalo_df is a dataframe of subhalo mass within snapshots.
# -----------------------------------------------------------------------------------------------------