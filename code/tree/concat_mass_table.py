import os
import argparse
from func import *
import pandas as pd
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default='Hydro-Arepo/MTNG-L500-4320-A/',  type=str)
args = parser.parse_args()

print('')
print(f'>>> Concatenate mass table <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

load_dir = f'result/DMhalo_mass_table/{args.sim}/'
fname_list = os.listdir(load_dir)
fname_list = [fname for fname in fname_list if fname.startswith('chunk')]

tot_df_list = []
for fname in tqdm(fname_list, desc='concatenate chunk df'):
    df = pd.read_csv(f'{load_dir}/{fname}', index_col=0)  # Assuming the first column is the index
    tot_df_list.append(df)
tot_df = pd.concat(tot_df_list, axis=1)


# Only keep necessary snapshots
final_index_values = ['snap_' + str(i) for i in [51, 69, 80, 94, 129, 151, 179, 214, 237, 264]]  
tot_df = tot_df.loc[final_index_values]

tot_df.to_csv(f'{load_dir}/mass_table.csv', index=True)