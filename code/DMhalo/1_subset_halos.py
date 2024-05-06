import illustris_python as il
import os
import argparse
import numpy as np

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--snap',default=78,type=int)
parser.add_argument('--mass_range',default=3,type=float)
args = parser.parse_args()

# Specify the snapshot
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'
snapNum = args.snap

# Select a subset of DM halos
Halos = il.groupcat.loadHalos(basePath, snapNum, fields='Group_M_Mean200')
subset_idx = np.where((Halos >= 10**args.mass_range) & 
                      (Halos < 10**(args.mass_range+0.5)))[0]
Ngroups_subset = subset_idx.shape[0]
print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.mass_range+10} ~ 10^{args.mass_range+10.5} MSun in at snap {snapNum}')

# Iterate over DM halos
for i, halo_idx in enumerate(subset_idx):
    if not os.path.exists(f'result/DM_halo_density_profiles/snap_{snapNum}/halo_{halo_idx}.npy'):
        os.system(f'python3 code/DMhalo/one_halo.py --snap {snapNum} --halo_idx {halo_idx}')
    else:
        print(f'At snap {snapNum}, DM halo local index {i+1}/{Ngroups_subset} already exists.')

print(f'All DM halos in the subset at snap {snapNum} are finished.')