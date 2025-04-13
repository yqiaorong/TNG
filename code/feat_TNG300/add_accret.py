import os
import numpy as np
import argparse
from tqdm import tqdm
import illustris_python as il
from func_accret import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='TNG300/sim_205_1250_Hydro', type=str)
parser.add_argument('--snapnum', default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add accretions to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Create a snapshot and redshift dict
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if 'DM' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'_DM/output/'
elif 'Hydro' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'/output/'
    
z_dict, a_dict = {}, {}
for snap in range(100):
    header = il.groupcat.loadHeader(basePath, snap)
    z_dict[snap] = header['Redshift']
    a_dict[snap] = header['Time']



# Load halo data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    print(data.keys())
    print(data['halo_M_Mean200'].shape)
    
    # Calculate one dynamical time from current z
    prev_snap = calc_tdyn(z_dict, data['z'])
    prev_masses = il.groupcat.loadHalos(basePath, prev_snap, fields='Group_M_Mean200') # Comoving mass!
    
    # Find accretion rates
    accretions = []
    for current_mass, subhalo_id in tqdm(zip(data['halo_M_Mean200'], data['FirstSub']), 
                                         total=len(data['FirstSub'])): # Physical mass!
        # Find the first subhalo history
        result = il.lhalotree.loadTree(basePath, args.snapnum, subhalo_id, fields=['SnapNum','SubhaloGrNr'], onlyMPB=True)
       
        if result is None:
            accretions.append([np.nan])
        else:
            # Only extract the result content at prev_snap
            target_idx = np.where(result['SnapNum'] == prev_snap)[0]
            if target_idx.size == 0:
                # print('No previous snapshot found! add nan')
                accretions.append([np.nan])
            else:
                prev_snap_SubhaloGrNr = result['SubhaloGrNr'][target_idx]
                # print(result['SnapNum'][target_idx], prev_snap_SubhaloGrNr)
                # Find the corresponding halo mass
                prev_mass = prev_masses[prev_snap_SubhaloGrNr] # Comoving mass!
                if prev_mass == 0:
                    # print('No previous mass found! add nan')
                    accretions.append([np.nan])
                else:
                    # Find the accretion rate
                    rate = np.log10(current_mass/(prev_mass/data['h'])) / np.log10(a_dict[args.snapnum]/a_dict[prev_snap]) # physical mass ratio! 
                    # print(rate)
                    accretions.append(rate)
    print(np.concatenate(accretions).shape)
    # Add accretion rates to data
    data['accretions'] = np.concatenate(accretions)
    # Save data
    np.save(halos_dir+fname, data)
    print('Saved! ')