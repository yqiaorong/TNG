import illustris_python as il
import numpy as np
import os
import argparse
import h5py
from tqdm import tqdm
from func import compt_density_profile

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--halo_idx',default=0,type=int)
args = parser.parse_args()

print('')
print(f'>>> DM halo density profile <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Specify the snapshot
basePath = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/L205n1250TNG/output'
snapNum = 99

# Load Halos from groupcat
group_fields = ['GroupCM', 'Group_M_Mean200', 'Group_R_Mean200']
Halos = il.groupcat.loadHalos(basePath, snapNum, fields=group_fields)

GroupCM = Halos['GroupCM']
Group_M_Mean200 = Halos['Group_M_Mean200']
Group_R_Mean200 = Halos['Group_R_Mean200']
del Halos

# Select DM halo
haloCM = GroupCM[args.halo_idx]
halo_M_Mean200 = Group_M_Mean200[args.halo_idx]
halo_R_Mean200 = Group_R_Mean200[args.halo_idx]
del GroupCM, Group_M_Mean200, Group_R_Mean200

# Set halo spatial boundary
xmin, xmax = haloCM[0] - 5 * halo_R_Mean200, haloCM[0] + 5 * halo_R_Mean200
ymin, ymax = haloCM[1] - 5 * halo_R_Mean200, haloCM[1] + 5 * halo_R_Mean200
zmin, zmax = haloCM[2] - 5 * halo_R_Mean200, haloCM[2] + 5 * halo_R_Mean200

# Load Halos coordinates from snapshot
# Coordinates = il.snapshot.loadSubset(basePath, snapNum, 'dm', ['Coordinates'], float32=True)

load_dir = os.path.join(basePath, 'snapdir_099')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.endswith('hdf5')]

rho_bins, r_bins = [], []
for idx, file in enumerate(load_list):
    
    # Load coordinates
    snap_path = os.path.join(load_dir, file)
    with h5py.File(snap_path, 'r') as f:

        coords = np.array(f['PartType1/Coordinates'])
        mask = (coords[:,0] >= xmin) & (coords[:,0] <= xmax) & (coords[:,1] >= ymin) & (coords[:,1] <= ymax) & (coords[:,2] >= zmin) & (coords[:,2] <= zmax)
        
        sub_coords = coords[mask]
        del coords, mask
    
    # Compute the density profile
    if sub_coords.shape[0] != 0:
        rho, r = compt_density_profile(sub_coords, haloCM, halo_R_Mean200)
        rho_bins.append(rho)
        r_bins = r
    del sub_coords
    
rho_bins = np.array(rho_bins)

# Get the final density profile
rho_bins = np.sum(rho_bins, axis=0)

# Save directory
save_dir = f'result/DM_halo_density_profiles/snap_{snapNum}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
# Save data
save_data = {'radius_profile': r_bins, 'density_profile': rho_bins, 
            'halo_R_Mean200': halo_R_Mean200, 'halo_M_Mean200': halo_M_Mean200}
np.save(os.path.join(save_dir, f'halo_{args.halo_idx}'), save_data)