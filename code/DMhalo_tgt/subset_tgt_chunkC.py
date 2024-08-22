import illustris_python as il
import h5py
import os
import argparse
import numpy as np
from tqdm import tqdm
from func import compt_density_profile
import time



### C func setup ###
import ctypes

# Load the lib in C
lib = ctypes.CDLL('code/DMhalo_tgt/dm_coords.so')

# Define the argument and return types for the C functions
lib.create_mask.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_bool), ctypes.c_int,
                            ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float]
lib.select_sub_coords.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_bool), ctypes.c_int]
lib.select_sub_coords.restype = ctypes.c_int



# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default=None,                     type=str)
parser.add_argument('--snapnum',      default=264,                      type=int)
parser.add_argument('--chunk_idx',    default=None,                      type=int)
parser.add_argument('--bin_start',    default=1,                        type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=5.5,                      type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',type=str)
args = parser.parse_args()

print('')
print('>>> DM halo density profiles per bin per chunk file <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Specify the snapshot
basePath = f'/virgotng/mpa/MTNG/{args.sim}'
snapnum = args.snapnum



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
    
# Convert to ctypes
Coordinates_c = Coordinates.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

num_coords = Coordinates.shape[0]


# Initialising the final data
global_idx = np.empty((Ngroups_subset))
halo_R_Mean200 = np.empty((Ngroups_subset))
halo_M_Mean200 = np.empty((Ngroups_subset))
densities = np.empty((Ngroups_subset, 85))
radii = np.empty((Ngroups_subset, 85))









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
    
    # Initial the final array
    mask = np.zeros(num_coords, dtype=np.bool_)
    mask_c = mask.ctypes.data_as(ctypes.POINTER(ctypes.c_bool))
    
    selected_coords = np.empty((num_coords, 3), dtype=np.float32)
    selected_coords_c = selected_coords.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    
    t1 = time.time()
    # Create mask
    lib.create_mask(Coordinates_c, mask_c, num_coords, xmin, xmax, ymin, ymax, zmin, zmax)
    t2 = time.time()
    
    # Select coordinates
    num_selected = lib.select_sub_coords(Coordinates_c, selected_coords_c, mask_c, num_coords)
    sub_coords = selected_coords[:num_selected]
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