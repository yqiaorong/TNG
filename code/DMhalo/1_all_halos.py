import illustris_python as il
import os
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--snap',default=78,type=int)
args = parser.parse_args()

# Specify the snapshot
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'
snapNum = args.snap

# Get the total number of DM halos
Header = il.groupcat.loadHeader(basePath, snapNum)
Ngroups_Total = Header['Ngroups_Total']

# Iterate over DM halos
for halo_idx in range(Ngroups_Total):
    if not os.path.exists(f'result/DM_halo_density_profiles/snap_{snapNum}/halo_{halo_idx}.npy'):
        os.system(f'python3 code/DMhalo/one_halo.py --snap {snapNum} --halo_idx {halo_idx}')
    else:
        print(f'At snap {snapNum}, DM halo local index {halo_idx+1}/{Ngroups_Total} already exists.')

print(f'All halos at snap {snapNum} are finished.')