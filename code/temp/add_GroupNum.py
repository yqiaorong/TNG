import os
import numpy as np
import argparse
import illustris_python as il

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='TNG300/sim_205_1250_Hydro', type=str)
parser.add_argument('--snapnum', default=None, type=int)
args = parser.parse_args()



# Select a subset of DM halos from groupcat
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if 'DM' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'_DM/output/'
elif 'Hydro' in args.sim:
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'/output/'
Group_M_Mean200 = il.groupcat.loadHalos(basePath, args.snapnum, fields='Group_M_Mean200')
halo_in_snap_index = np.where((Group_M_Mean200 >= 10**1) & (Group_M_Mean200 < 10**4.5))[0] # Comoving mass
print(halo_in_snap_index.shape)


# Load halo data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    print(data['halo_M_Mean200'].shape)
    
    # Check if mass matches
    print(data['halo_M_Mean200'].shape)
    print(Group_M_Mean200[halo_in_snap_index].shape)
    if np.array_equal(sorted(data['halo_M_Mean200']), sorted(Group_M_Mean200[halo_in_snap_index]/data['h'])):
        print('Mass matches!')
    else:
        print('Mass does not match!')
        
    # Correct sorting
    sorted_halo_idx = np.argsort(data['halo_M_Mean200'])
    sorted_group_idx = np.argsort(Group_M_Mean200[halo_in_snap_index])
    
    correct_groupNum = halo_in_snap_index[sorted_group_idx[np.argsort(sorted_halo_idx)]]
    
    # Double check
    if np.array_equal(Group_M_Mean200[correct_groupNum]/data['h'], data['halo_M_Mean200']):
        print('Halo index correct!')
    else:
        print('Halo index incorrect!')
        os._exit(0)
        
    # Add halo index to data
    data['GroupNum'] = correct_groupNum
    data['FirstSub'] = il.groupcat.loadHalos(basePath, args.snapnum, fields='GroupFirstSub')[correct_groupNum]
    
    print(data['halo_M_Mean200'].shape, data['GroupNum'].shape, data['FirstSub'].shape)
    print(data.keys())
    # Save data
    np.save(halos_dir+fname, data)