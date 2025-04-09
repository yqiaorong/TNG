"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import os
import numpy as np
from DMhalo_TNG300.func import *
from tqdm import tqdm
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='MTNG/Hydro-Arepo/MTNG-L500-4320-A/', type=str)
parser.add_argument('--snapnum', default=None, type=int)
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
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
    densities      = data['densities']      # [Msun / (kpc)^3]
    radial_bins    = data['radial_bins']    # [kpc]
    # Compute the NSW profile
    all_rho0, all_Rs = [], []
    # idx = 0
    for r, rho, R200 in tqdm(zip(radial_bins, densities, halo_R_Mean200)):
        # Crop region r < R200
        rho = rho[r < R200]
        r = r[r < R200]
        # Remove zeros in densities
        non_zero_indices = np.where(rho > 0)
        r = r[non_zero_indices]
        rho = rho[non_zero_indices]        
        # Fit the log NSW profile
        log_rho = np.log10(rho)
        
        # Filter the density profile
        if args.sim.startswith('MTNG'):
            rsoft = 100 # [kpc]
            start_idx = np.where(r > rsoft)[0][0]
            r, log_rho= r[start_idx:], log_rho[start_idx:]
        
        popt, _ = fit_log_NFW_profile(r, log_rho, R200, 
                                    #   idx
                                      )

        rho_0, R_s = popt
        all_rho0.append(rho_0)
        all_Rs.append(R_s)
        # idx += 1
        
    all_rho0 = np.array(all_rho0)
    all_Rs = np.array(all_Rs)
    # Compute the concentrations
    conc = halo_R_Mean200 / all_Rs

    # Append data
    data['NFWrho0'] = all_rho0
    data['NFWRs']   = all_Rs
    data['NFWconc'] = conc
    
    # Save the data
    print(data.keys())
    np.save(halos_dir+fname, data)