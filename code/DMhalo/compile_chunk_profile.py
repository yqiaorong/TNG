from tqdm import tqdm
import numpy as np
from func import *
import illustris_python as il
import argparse
import os
import h5py

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',       default='_DM',type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--bin_start',default=1,    type=float)
parser.add_argument('--bin_end',  default=None, type=float)
args = parser.parse_args()

print('')
print(f'>>> BCompile chunk files of DM density profiles <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG'%(boxsize,res)+f'{args.DM}/output'



# Make mass cuts
bin_start = args.bin_start
bin_end = args.bin_end
print(f'The current mass range: 10^{bin_start+10} ~ 10^{bin_end+10} MSun/h')


    
snap = args.snapnum
# Load redshift values
with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    z = 1 / scale_factor - 1
    h = header['HubbleParam']
    
    
### Load halos ###

# Halos data root dir
halos_dir = f'result/{root_dir}/sim_{boxsize}_{res}{args.DM}/snap_{snap}/densities'

# First round of rough selection of halos based on M200
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snap, fields='Group_M_Mean200')
subset_idx = np.where((Group_M_Mean200 >= 10**bin_start) & (Group_M_Mean200 < 10**bin_end))[0]
halos_list = [f'halo_{idx}.npy' for idx in subset_idx]
print(f'The total halo numbers: {len(halos_list)}')
del Group_M_Mean200

if len(halos_list) == 0:
    pass 
else:
    total_R200, total_M200, total_rho, total_r = [], [], [], []
    ### Compile density profiles 
    for ihalo, halo in enumerate(tqdm(halos_list)):
        data = np.load(f'{halos_dir}/{halo}', allow_pickle=True).item()
        
        total_R200.append(data['halo_R_Mean200'])
        total_M200.append(data['halo_M_Mean200'])
        total_rho.append(data['densities'])
        total_r.append(data['radial_bins'])

    total_R200 = np.array(total_R200)
    total_M200 = np.array(total_M200)
    total_rho = np.array(total_rho)
    total_r = np.array(total_r)

    print(total_rho.shape, total_r.shape, total_R200.shape, total_M200.shape, subset_idx.shape)



    ### Save the compiled data
    save_dict = {'halo_R_Mean200': total_R200, # [ckpc/h]
                'halo_M_Mean200': total_M200, # [10^10 Msun/h]
                'densities': total_rho, # [(Msun/h)/(ckpc/h)^3]
                'radial_bins': total_r, # [ckpc/h]
                'h': h, 'scale_factor': scale_factor, 'z': z}  

    save_dir = f'result/DMhalo_density_profiles/TNG300/sim_{boxsize}_{res}{args.DM}/snap_{snap}/final_densities/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    np.save(f'{save_dir}/bin-{int(bin_start*10)}-{int(bin_end*10)}', save_dict) 
    print('saved')