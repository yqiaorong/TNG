import argparse
import numpy as np
import os
import h5py
import illustris_python as il
from matplotlib import pyplot as plt

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',   default=205,type=int)
parser.add_argument('--res',       default=1250,type=int)
parser.add_argument('--root_dir',  default=None, type=str)
args = parser.parse_args()

print('')
print(f'>>> Sort splashback radius data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# BasePath
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
boxsize = args.boxsize
res = args.res
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Load dir
load_dir = f'result/Stacked_{args.root_dir}/sim_{boxsize}_{res}'

# Walk through the directory
file_list = []
for dirpath, subfolder, fnames in os.walk(load_dir):
    for fname in fnames:
        if fname.endswith('.npy'):
            file_path = os.path.join(dirpath, fname.replace('-', ''))
            file_list.append(file_path)

# Sort index    
snap_list = sorted(os.listdir(load_dir))
file_list = sorted(file_list)



# Save data path
data_path = f'data/{args.root_dir}/sim_{boxsize}_{res}'
if not os.path.exists(data_path):
    os.makedirs(data_path)
    
# Create hdf5 file
with h5py.File(os.path.join(data_path,'DMhalo_profiles.hdf5'), 'w') as data_f:
    
    # Iterate over snapshots
    for snap in snap_list:
        print(snap)
        # Create group 
        group = data_f.create_group(snap)
        
        # Select files within this snapshot
        subfile_list = [fname for fname in file_list if snap in fname]
        subfile_list = sorted(subfile_list)
        
        # Load redshift values
        with h5py.File(il.snapshot.snapPath(basePath, snap[-2:]), 'r') as f:
            header = dict(f['Header'].attrs.items())
            scale_factor = header['Time']
            z = 1 / scale_factor - 1
            h = header['HubbleParam']
        
        # Iterate over subfiles  
        R200_one_snap, rf_one_snap, pdf_one_snap, sf_one_snap = [], [], [], []
        for file in subfile_list:
            print(file)
            # Load data
            data = np.load(file, allow_pickle=True).item()
            
            # R200
            R200_median = data['R200_median']
            R200_one_snap.append(R200_median) # [ckpc/h]
            # Radius
            rf_one_snap.append(data['profile_radius_fit'] * R200_median * scale_factor / h) # [kpc]
            # Fitted density profiles
            pdf_one_snap.append(data['profile_densities_fit'])
            # Fitted gradients
            sf_one_snap.append(data['slopes_fit'])
        
        # Add data to the group
        group.create_dataset('z', data=z)
        group.create_dataset('R200', data=np.array(R200_one_snap))
        group.create_dataset('radius_fit', data=np.array(rf_one_snap))
        group.create_dataset('densities_fit', data=np.array(pdf_one_snap))
        group.create_dataset('slopes_fit', data=np.array(sf_one_snap))  
        print('')  