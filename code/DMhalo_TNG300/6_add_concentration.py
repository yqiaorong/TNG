"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import os
import numpy as np
from func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default='TNG300/sim_205_1250_Hydro/', type=str)
parser.add_argument('--snapnum',  default=13, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add concentrations to halo data <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Load halos data
# -----------------------------------------------------------------------------------------
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    print(fname)
    data = np.load(os.path.join(halos_dir, fname), allow_pickle=True).item()
    # z            = data['z']
    # scale_factor = data['scale_factor']
    # h            = data['h']
    # rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
    
    halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
    densities      = data['densities']      # [Msun / (kpc)^3]
    radial_bins    = data['radial_bins']    # [kpc]
    del data
    
    # Compute the NSW profile
    for r, rho, R200 in zip(radial_bins[:10], densities[:10], halo_R_Mean200[:10]):
        popt, _ = fit_NSW_profile(r, rho, R200)
