import illustris_python as il
import os
import argparse
import numpy as np

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',    default=1250,type=int)
parser.add_argument('--DM',    default='_DM',type=str)
parser.add_argument('--bin_start',default=1, type=float)
parser.add_argument('--bin_end',default=4.5, type=float)
args = parser.parse_args()

print('')
print(f'>>> Check halo number <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG'%(args.boxsize,args.res)+f'{args.DM}/output/'

# Select a subset of DM halos
snaps = [# 8, 13, 17, 21, 25, 
         33, 40, 50, 67, 78, 99]
for s in snaps:
    Halos = il.groupcat.loadHalos(basePath, s, fields='Group_M_Mean200')
    subset_idx = np.where((Halos >= 10**args.bin_start) & 
                        (Halos < 10**args.bin_end))[0]
    Ngroups_subset = subset_idx.shape[0]
    print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.bin_start+10} ~ 10^{args.bin_end+10} MSun/h in at snap {s}')
    
    halos_set = set(os.listdir(f'result/DMhalo_density_profiles_old/sim_{args.boxsize}_{args.res}{args.DM}/snap_{s}/densities'))
    subset_idx_set = {f'halo_{idx}.npy' for idx in subset_idx}
    num_missing = len(subset_idx_set.difference(halos_set))
    print(f'{num_missing} halos missing.')