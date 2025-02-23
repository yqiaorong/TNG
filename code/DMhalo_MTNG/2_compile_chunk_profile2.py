import os
import argparse
import numpy as np
from tqdm import tqdm

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',          default='DM-Arepo/MTNG-L500-4320-A',type=str)
parser.add_argument('--snapnum',      default=264,                        type=int)
parser.add_argument('--bin_start',    default=2,                          type=float) # [10^{10+x} Msun/h]
parser.add_argument('--bin_end',      default=3,                          type=float) # [10^{10+x} Msun/h]
parser.add_argument('--save_root_dir',default='DMhalo_density_profiles',  type=str)
args = parser.parse_args()

print('')
print('>>> Compile DM halo density profiles per mass cut <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



### Load intermediate density profiles ###
parent_dir = f'/nfs/mvogelsblab001/Users/s_qyu/{args.save_root_dir}/{args.sim}/snap_{args.snapnum}/'
load_dir = f'{parent_dir}/intermediate_densities/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}/'
load_list = os.listdir(load_dir)
print(len(load_list))

if len(load_list) != 640:
    os._exit(0)
else:
    # Get halo numbers in the mass cut
    data = np.load(f'{load_dir}/{load_list[0]}', allow_pickle=True).item()
    num_halos = data['halo_R_Mean200'].shape[0]

    # Load data
    densities = np.zeros((num_halos, 85))
    for ifname, fname in enumerate(tqdm(load_list, desc='chunk files')):
        data = np.load(f'{load_dir}/{fname}', allow_pickle=True).item()
        if ifname == 0: 
            halo_R_Mean200 = data['halo_R_Mean200'] # [ckpc / h]
            halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 Msun / h]
            radii = data['radial_bins']             # [ckpc / h]
            h = data['h']
            scale_factor = data['scale_factor']
            z = 1 / scale_factor - 1
        densities += data['densities']
        del data
    print(densities.shape, radii.shape, halo_R_Mean200.shape, halo_M_Mean200.shape)



    ### Save compiled chunk file ###
    save_dict = {'halo_R_Mean200': halo_R_Mean200, # [ckpc/h]
                'halo_M_Mean200': halo_M_Mean200,  # [10^10 Msun/h]
                'densities': densities,            # [(Msun/h)/(ckpc/h)^3]
                'radial_bins': radii,              # [ckpc/h]
                'h': h, 'scale_factor': scale_factor, 'z': z}   

    # Save 
    save_dir = f'{parent_dir}/final_densities'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    np.save(f'{save_dir}/bin-{int(args.bin_start*10)}-{int(args.bin_end*10)}', save_dict) 