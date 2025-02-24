import os
import numpy as np
import illustris_python as il
from tqdm import tqdm
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim_type', default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Mass table 2 <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

sim_type, snapnum = args.sim_type, args.snapnum



# Load halo density profiles
# ------------------------------------------------------------------------------
halos_dir = f'result/DMhalo_density_profiles_phys/MTNG/{sim_type}-Arepo/MTNG-L500-4320-A/snap_{snapnum}/final_densities/'
halos_list = os.listdir(halos_dir)
# sort the list
halos_list = sorted(halos_list)
print(halos_list)
# Get the bin starts and ends
bin_starts = [float(halo.split('-')[1])/10 for halo in halos_list]
bin_ends = [float(halo.split('-')[2].split('.')[0])/10 for halo in halos_list]
print(bin_starts)
print(bin_ends)  

# Load all halo masses
# ------------------------------------------------------------------------------
sim = f'{sim_type}-Arepo/MTNG-L500-4320-A'
basePath = f'/virgotng/mpa/MTNG/{sim}/output/'
snap_all_mass = il.groupcat.loadHalos(basePath, snapnum, fields='Group_M_Mean200')
print(snap_all_mass.shape)

# Load FPG mass, idx, and linking subhalo mass
# ------------------------------------------------------------------------------
table_dir = f'result/DMhalo_mass_table_new/{sim_type}-Arepo/MTNG-L500-4320-A/'
# This is physical mass!
FPGrMass = np.load(table_dir+f'snap_{snapnum}_FPGrMass.npy')
print(FPGrMass.shape)
FPGr = np.load(table_dir+f'snap_{snapnum}_FPGr.npy')
print(FPGr.shape)
SubMass = np.load(table_dir+f'snap_{snapnum}_SubMass.npy')
print(SubMass.shape)

# Load accretion rates
# ------------------------------------------------------------------------------
accret_data = np.load(f'result/DMhalo_mass_table_new/{sim_type}-Arepo/MTNG-L500-4320-A/accretion_rates.npy',
                          allow_pickle=True).item()
snap_idx = np.where(accret_data['snaps'] == np.array(snapnum))[0][0]
accret_rates = accret_data['accretions'][:, snap_idx]
print(accret_rates.shape)
del accret_data

# Add accretions to halo density profiles
# ------------------------------------------------------------------------------
halo_M, halo_R, bins, densities, accretions = [], [], [], [], []
for fname, start, end in zip(halos_list, bin_starts, bin_ends):
    print(fname)
    
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    halo_M.append(data['halo_M_Mean200']) # This is physical mass!
    halo_R.append(data['halo_R_Mean200'])
    bins.append(data['radial_bins'])
    densities.append(data['densities'])
    print(data['halo_M_Mean200'].shape)
    
    # Get the halo snap indices
    subset_idx = np.where((snap_all_mass >= 10**start) & (snap_all_mass < 10**end))[0]
    print(subset_idx.shape)
    
    # Add the accretion rates
    for idx in tqdm(subset_idx):
        # Check the masses
        idx_in_FPGr = np.where(FPGr == idx)[0]
        # print(FPGrMass[idx_in_FPGr].shape, snap_all_mass[idx])
        # print(np.all(FPGrMass[idx_in_FPGr] == FPGrMass[idx_in_FPGr][0]), FPGrMass[idx_in_FPGr][0])
        
        if FPGrMass[idx_in_FPGr].shape[0] == 0:
            accretions.append([np.nan])
        else:
            if np.all(FPGrMass[idx_in_FPGr] == FPGrMass[idx_in_FPGr][0]) and FPGrMass[idx_in_FPGr][0] == snap_all_mass[idx]:
                if idx_in_FPGr.shape[0] == 1:
                    right_accret = accret_rates[idx_in_FPGr[0]]
                else:
                    # Select the one with the highest subhalo mass
                    idx_in_FPGr = idx_in_FPGr[np.argmax(SubMass[idx_in_FPGr])]
                    right_accret = accret_rates[idx_in_FPGr]
                accretions.append([right_accret])
            else:
                exit()

# Concatenate data
# ------------------------------------------------------------------------------
halo_M     = np.concatenate(halo_M)
halo_R     = np.concatenate(halo_R)
bins       = np.concatenate(bins)
densities  = np.concatenate(densities)
accretions = np.concatenate(accretions)
print(halo_M.shape, halo_R.shape, bins.shape, densities.shape, accretions.shape)
    
# Saved dict
save_dict = {'h': data['h'], 
            'scale_factor': data['scale_factor'], 
            'z': data['z'], 
            'rho_c': data['rho_c'],
            'halo_M_Mean200': halo_M,
            'halo_R_Mean200': halo_R,
            'radial_bins':    bins,
            'densities':      densities,
            'accretion_rate': accretions}

# Create the save path
save_dir = f'result/DMhalo_density_profiles/MTNG/{sim_type}-Arepo/MTNG-L500-4320-A/snap_{snapnum}/final_densities/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
np.save(os.path.join(save_dir, f'bin-{int(bin_starts[0]*10)}-{int(bin_ends[-1]*10)}'), save_dict)