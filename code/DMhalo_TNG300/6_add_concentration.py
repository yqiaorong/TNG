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


halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    # Load data
    data = np.load(os.path.join(halos_dir, fname), allow_pickle=True).item()
    halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
    densities      = data['densities']      # [Msun / (kpc)^3]
    radial_bins    = data['radial_bins']    # [kpc]
    # Compute the NSW profile
    all_rho0, all_Rs = [], []
    for r, rho, R200 in zip(radial_bins, densities, halo_R_Mean200):
        # Crop region r < R200
        rho = rho[r < R200]
        r = r[r < R200]
        # Remove zeros in densities
        non_zero_indices = np.where(rho > 0)
        r = r[non_zero_indices]
        rho = rho[non_zero_indices]
        # Fit the NSW profile
        popt, _ = fit_NSW_profile(r, rho, R200)
        rho_0, R_s = popt
        all_rho0.append(rho_0)
        all_Rs.append(R_s)
    all_rho0 = np.array(all_rho0)
    all_Rs = np.array(all_Rs)
    # Compute the concentrations
    conc = all_Rs / halo_R_Mean200
    # Append data
    data['NSW_rho0'] = all_rho0
    data['NSW_Rs'] = all_Rs
    data['NSW_conc'] = conc
    # Save the data
    print(data.keys())
    # np.save(os.path.join(halos_dir, fname), data)