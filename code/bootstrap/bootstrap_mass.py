import os
import numpy as np
from boots_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1024, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per mass cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load halos data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fname = os.listdir(halos_dir)[0]
print(halos_fname)
data = np.load(os.path.join(halos_dir, halos_fname), allow_pickle=True).item()

z            = data['z']
scale_factor = data['scale_factor']
h            = data['h']
rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 Msun]
densities      = data['densities']      # [Msun / (kpc)^3]
radial_bins    = data['radial_bins']    # [kpc]
del data

print(f'total number of halos: {halo_M_Mean200.shape[0]}')


# Bootstrap setup
final_results = bootstrap_all_features(args, 
                                [z, h, rho_c], 
                                [radial_bins, densities, halo_M_Mean200, halo_R_Mean200])     
    
    
# Save the results
save_stats_dir = f'result/bootstrap_stats/with_mass/{args.sim}/Nboots_{args.Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)

np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', final_results)