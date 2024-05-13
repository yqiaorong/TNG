import illustris_python as il
import os
import argparse
import numpy as np

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
parser.add_argument('--snapnum',default=99,type=int)
parser.add_argument('--mass_range',default=3,type=float)
args = parser.parse_args()

# Specify the snapshot
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(args.boxsize,args.res)
snapnum = args.snapnum

# Select a subset of DM halos
Halos = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')
subset_idx = np.where((Halos >= 10**args.mass_range) & 
                      (Halos < 10**4))[0]
Ngroups_subset = subset_idx.shape[0]
print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.mass_range+10} ~ 10^14 MSun in at snap {snapnum}')

# Iterate over DM halos
missing_number = 0
for i, halo_idx in enumerate(subset_idx):
    if not os.path.exists(f'result/DMhalo_density_profiles/sim_{args.boxsize}_{args.res}/snap_{snapnum}/densities/halo_{halo_idx}.npy'):
        missing_number += 1
print(f'{missing_number} halos missing.')