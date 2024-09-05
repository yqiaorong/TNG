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
print(f'>>> Coordinates table <<<')
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



index_values = ['snap_' + str(i) for i in range(2, 99)]
final_index_values = ['snap_' + str(i) for i in [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]]  

# Save directory
save_dir = f'result/DMhalo_idx_table/sim_{boxsize}_{res}{args.DM}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Concatenate all dataframes
df_list = os.listdir(save_dir)[:2] # Just for test so didn't use all the halos to save time
tot_df_list = []
for file in tqdm(df_list, desc='concatenate chunk df'):
    df = pd.read_csv(f'{save_dir}/{file}', index_col=0)  # Assuming the first column is the index
    tot_df_list.append(df)
tot_df = pd.concat(tot_df_list, axis=1)

df_y = tot_df.copy()
df_z = tot_df.copy()

# Save the dataframe
save_mass_dir = f'result/DMhalo_coords_table/sim_{boxsize}_{res}{args.DM}'
if not os.path.exists(save_mass_dir):
    os.makedirs(save_mass_dir)

# Now change the entries of these global index to their masses
for idf, df in enumerate([tot_df, df_y, df_z]):
    df = get_field_values_of_lifeline(basePath, df, group_field='GroupPos', coords_idx=0)
    df.to_csv(f'{save_mass_dir}/coords_{idf}_table.csv', index=final_index_values)