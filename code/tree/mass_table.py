import os
import h5py
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
parser.add_argument('--DM',     default='',    type=str)
parser.add_argument('--bin_start',default=1,type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',default=4.5,type=float) # [10^{10+x} Msun/h]
args = parser.parse_args()

print('')
print(f'>>> Mass table <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG'%(args.boxsize,args.res)+f'{args.DM}/'
snapnum = 99 # The last snapshot, should not change
boxsize = args.boxsize
res     = args.res


# Select a subset of DM halos from groupcat
group_fields = ['Group_M_Mean200', 'GroupFirstSub']
Halos = il.groupcat.loadHalos(basePath+'output', snapnum, fields=group_fields)

Group_M_Mean200 = Halos['Group_M_Mean200'] # [10^10 MSun/h]
GroupFirstSub = Halos['GroupFirstSub']

# Using physical mass to select subset
halo_global_idx = np.where((Group_M_Mean200 >= 10**args.bin_start) & 
                           (Group_M_Mean200 < 10**args.bin_end))[0]
num_halos = halo_global_idx.shape[0]
print(f'The number of halos in the given mass range: {num_halos}')

index_values = ['snap_' + str(i) for i in range(2, 99)]
final_index_values = ['snap_' + str(i) for i in [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]]  

# Save directory
save_dir = f'result/DMhalo_idx_table/sim_{boxsize}_{res}{args.DM}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

chunk_size = 100
chunk_indices = range(0, num_halos, chunk_size)
for chunk_idx in chunk_indices:
    if not os.path.exists(f'{save_dir}/chunk_{chunk_idx}.csv'):
        ### Creat halo (first subhalo) lifeline dataframe ###
        df = pd.DataFrame(index=index_values)

        # Load treeX 
        for idx in tqdm(halo_global_idx[chunk_idx: min(num_halos, chunk_idx+chunk_size)], 
                        desc=f'make index df: chunk {chunk_idx}+'):
            subhalo_global_idx = GroupFirstSub[idx]

            # Get the dataframe of global index
            lifeline(basePath, subhalo_global_idx, df)

        # Only keep necessary snapshots
        df = df.loc[final_index_values]

        # Save the dataframe
        df.to_csv(f'result/DMhalo_idx_table/sim_{boxsize}_{res}{args.DM}/chunk_{chunk_idx}.csv', index=True)

# Concatenate all dataframes
df_list = os.listdir(save_dir)
tot_df_list = []
for file in tqdm(df_list, desc='concatenate chunk df'):
    df = pd.read_csv(f'{save_dir}/{file}', index_col=0)  # Assuming the first column is the index
    tot_df_list.append(df)
tot_df = pd.concat(tot_df_list, axis=1)

# Now change the entries of these global index to their masses
tot_df = get_field_values_of_lifeline(basePath, tot_df, group_field='Group_M_Mean200')

# Save the dataframe
save_mass_dir = f'result/DMhalo_mass_table/sim_{boxsize}_{res}{args.DM}'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)
tot_df.to_csv(f'{save_mass_dir}/mass_table.csv', index=final_index_values)