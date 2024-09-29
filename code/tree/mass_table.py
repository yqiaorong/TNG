import os
import argparse
from func import *
import pandas as pd

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',    default='Hydro-Arepo/MTNG-L500-4320-A/',  type=str)
parser.add_argument('--chunk_start_idx', default=None,  type=int)
args = parser.parse_args()

print('')
print(f'>>> Mass table <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



basePath = f'/virgotng/mpa/MTNG/{args.sim}/output/'

# Save directory
save_dir = f'result/DMhalo_mass_table/{args.sim}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)



# Iterate over chunk files
for chunk_idx in range(args.chunk_start_idx, 1):
    print('chunk idx: ', chunk_idx)
    if not os.path.exists(f'{save_dir}/chunk_{chunk_idx}.csv'):

        # Get the dataframe of Group Mass
        df = Lifeline2(basePath+'treedata/', chunk_idx)
        if df.columns.shape[0] == 0:
            pass
        else:
            # Save the chunk dataframe
            df.to_csv(f'{save_dir}/chunk_{chunk_idx}.csv', index=True)