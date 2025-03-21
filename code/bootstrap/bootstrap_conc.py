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
print(f'>>> Bootstrap splashback features per conc cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load halos data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    data = np.load(os.path.join(halos_dir, fname), allow_pickle=True).item()

    z            = data['z']
    scale_factor = data['scale_factor']
    h            = data['h']
    rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
    halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
    NFW_conc       = data['NFW_conc']
    densities      = data['densities']      # [Msun / (kpc)^3]
    radial_bins    = data['radial_bins']    # [kpc]
    del data

# Bootstrap setup
final_results = bootstrap_all_features(args, 
                                [z, h, rho_c], 
                                [radial_bins, densities, NFW_conc, halo_R_Mean200],
                                'conc', bin_start=0, bin_end=5, bin_width=0.2)     
    
    
# Save the results
save_stats_dir = f'result/bootstrap_stats/with_conc/{args.sim}/Nboots_{args.Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)

np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', final_results)