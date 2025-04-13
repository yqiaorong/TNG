from tqdm import tqdm
import numpy as np
import illustris_python as il
import argparse
import os
import h5py

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',       default='Hydro',type=str)
parser.add_argument('--snapnum',  default=None, type=int)
args = parser.parse_args()

print('')
print('>>> Compile chunk files of DM halo density profiles <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



root_dir = 'DMhalo_density_profiles_raw'
boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
if args.DM == 'DM':
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'_DM/output/'
elif args.DM == 'Hydro':
    basePath = data_path + 'L%dn%dTNG'%(205, 1250)+'/output/'



# Make mass cuts
bin_start = 1
bin_end = 5
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
halos_dir = f'result/{root_dir}/sim_{boxsize}_{res}_{args.DM}/snap_{snap}/densities'
halos_list = os.listdir(halos_dir)
print(f'The total halo numbers: {len(halos_list)}')

if len(halos_list) == 0:
    pass 
else:
    total_R200, total_M200, total_rho, total_r = [], [], [], []
    total_GrNr, total_FirstSub = [], []
    ### Compile density profiles 
    for halo in tqdm(halos_list):
        data = np.load(f'{halos_dir}/{halo}', allow_pickle=True).item()
        
        total_R200.append(data['halo_R_Mean200'])
        total_M200.append(data['halo_M_Mean200'])
        total_rho.append(data['densities'])
        total_r.append(data['radial_bins'])
        total_GrNr.append(data['GroupNum'])
        total_FirstSub.append(data['GroupFirstSub'])

    total_R200 = np.array(total_R200)
    total_M200 = np.array(total_M200)
    total_rho = np.array(total_rho)
    total_r = np.array(total_r)
    total_GrNr = np.array(total_GrNr)
    total_FirstSub = np.array(total_FirstSub)

    print(total_rho.shape, total_r.shape, total_R200.shape, total_M200.shape, total_GrNr.shape, total_FirstSub.shape)



    ### Save the compiled data
    save_dict = {'halo_R_Mean200': total_R200, # [ckpc/h]
                'halo_M_Mean200':  total_M200, # [10^10 Msun/h]
                'densities':       total_rho,  # [(Msun/h)/(ckpc/h)^3]
                'radial_bins':     total_r,    # [ckpc/h]
                'GroupNum':        total_GrNr,
                'FirstSub':        total_FirstSub,
                'h': h, 'scale_factor': scale_factor, 'z': z}  
    print(save_dict.keys())

    save_dir = f'result/DMhalo_density_profiles_raw2/TNG300/sim_{boxsize}_{res}_{args.DM}/snap_{snap}/final_densities/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    np.save(f'{save_dir}/bin-{int(bin_start*10)}-{int(bin_end*10)}', save_dict) 
    print('saved')