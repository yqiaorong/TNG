import illustris_python as il
import h5py
import os
import argparse
import numpy as np
from tqdm import tqdm
from func import compt_density_profile
import time
from numba import njit,prange

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default=None,                     type=str)
parser.add_argument('--snapnum',      default=264,                      type=int)
parser.add_argument('--chunk_idx',    default=None,                     type=int)
parser.add_argument('--bin_start',    default=3,                        type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=3.5,                      type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

print('')
print('>>> DM halo density profiles per bin per chunk file <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
snapnum = args.snapnum
basePath = f'/nfs/mvogelsblab001/Users/s_qyu/MTNG/{args.sim}/'



# Save root dir
save_root_dir = args.save_root_dir
    
    

# Select a subset of DM halos from groupcat
group_fields = ['GroupPos', 'Group_M_Mean200', 'Group_R_Mean200']
Halos = il.groupcat.loadHalos(basePath, snapnum, fields=group_fields)

GroupPos        = Halos['GroupPos']        # [ckpc/h]
Group_M_Mean200 = Halos['Group_M_Mean200'] # [10^10 MSun/h]
Group_R_Mean200 = Halos['Group_R_Mean200'] # [ckpc/h]

# Using physical mass to select subset
subset_idx = np.where((Group_M_Mean200 >= 10**args.bin_start) & (Group_M_Mean200 < 10**args.bin_end))[0]
Ngroups_subset = subset_idx.shape[0]
print(f'In total, {Ngroups_subset} DM halos with mass 10^{args.bin_start+10} ~ 10^{args.bin_end+10} MSun in at snap {snapnum}')



# Load params
with h5py.File(il.snapshot.snapPath(basePath, snapnum), 'r') as f:
    header = dict(f['Header'].attrs.items())
    Parameters = dict(f['Parameters'].attrs.items())

    BoxSize = Parameters['BoxSize'] # [cMpc / h]
    scale_factor = header['Time']
    h = Parameters['HubbleParam']
    DMmass = header['MassTable'][1] * 10**10 # [MSun / h]
    


# Load coordinates ONLY AT ONCE
load_dir = os.path.join(basePath, f'snapdir_{snapnum:03d}')
load_list = os.listdir(load_dir)
load_list = [fname for fname in load_list if fname.startswith(f'snapshot_{snapnum:03d}') and fname.endswith('hdf5')]

# Load coordinates
snap_path = os.path.join(load_dir, load_list[args.chunk_idx])
with h5py.File(snap_path, 'r') as f:
    Coordinates = np.array(f['PartType1/Coordinates'], dtype='float32') # [ckpc/h]
    print(Coordinates.shape)



# Initialising the final data
global_idx     = np.empty((Ngroups_subset))
halo_R_Mean200 = np.empty((Ngroups_subset))
halo_M_Mean200 = np.empty((Ngroups_subset))
densities      = np.empty((Ngroups_subset, 85))
radii          = np.empty((Ngroups_subset, 85))


### SPEED UP!!! ###

@njit(parallel=True)
def create_mask(Coordinates, xmin, xmax, ymin, ymax, zmin, zmax):
    num_coords = Coordinates.shape[0]
    mask = np.zeros(num_coords, dtype=np.bool_)
    for i in prange(num_coords):
        x, y, z = Coordinates[i]
        if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
            mask[i] = True
    return mask

@njit(parallel=True)
def select_sub_coords(Coordinates, mask):
    num_coords = Coordinates.shape[0]
    selected_coords = np.empty((np.sum(mask), Coordinates.shape[1]), dtype=Coordinates.dtype)
    
    idx = 0
    for i in prange(num_coords):
        if mask[i]:
            selected_coords[idx] = Coordinates[i]
            idx += 1
    
    return selected_coords



# Iterate over DM halos
ti = time.time()
for i, idx in enumerate(tqdm(subset_idx, desc=f'chunk {args.chunk_idx}')):
    
    # Initial values 
    x, y, z = np.round(GroupPos[idx, 0].item(), 1), np.round(GroupPos[idx, 1].item(), 1), np.round(GroupPos[idx, 2].item(), 1)
    R = np.round(Group_R_Mean200[idx].item(), 3)

    # Set halo spatial boundary
    edge = 5 # as a multiple of hal0_R_Mean200
    xmin, xmax = x - edge * R, x + edge * R
    ymin, ymax = y - edge * R, y + edge * R
    zmin, zmax = z - edge * R, z + edge * R
    
    # Select region of coordinates
    t1 = time.time()
    mask = create_mask(Coordinates, xmin, xmax, ymin, ymax, zmin, zmax)
    t2 = time.time()
    sub_coords = select_sub_coords(Coordinates, mask)
    t3 = time.time()
    print(f'{t2-t1:.4f}, {t3-t2:.4f}')
    
    # Compute the density profile
    rhos, radial_bins = compt_density_profile(sub_coords, [x,y,z], R, BoxSize*1000)
    rhos = rhos * DMmass
    
    # Save the density profile
    global_idx[i] = idx
    halo_R_Mean200[i] = Group_R_Mean200[idx]
    halo_M_Mean200[i] = Group_M_Mean200[idx]
    densities[i] = rhos
    radii[i] = radial_bins
tf = time.time()
print(f'The total time: {tf-ti:.4f}')

print(global_idx.shape, halo_R_Mean200.shape, halo_M_Mean200.shape)
print(densities.shape, radii.shape)
save_dict = {'halo_R_Mean200': halo_R_Mean200, # [ckpc / h]
             'halo_M_Mean200': halo_M_Mean200, # [10^10 Msun / h]
             'h': h, 'scale_factor': scale_factor, 'DMmass': DMmass, # [MSun / h]
             'densities': densities,            # [(Msun/h)/(ckpc/h)^3]
             'radial_bins': radii               # [ckpc/h]
            }     



# Save directory
save_dir = 'result/'+args.save_root_dir+f'/{args.sim}/snap_{snapnum}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    
save_data_dir = f'{save_dir}/intermediate_densities/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}'
if not os.path.exists(save_data_dir):
    os.makedirs(save_data_dir)
    
# Save
np.save(os.path.join(save_data_dir, f'chunk-{args.chunk_idx}'), save_dict) 
print(f'chunk {args.chunk_idx}: saved')