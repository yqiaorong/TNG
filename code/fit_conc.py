"""All the data computed in this script is saved in result/DMhalo_density_profiles/"""

import os
import numpy as np
from tqdm import tqdm
import argparse
from bootstrap.func import *
from splashback import splashback_properties 

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default='MTNG/Hydro-Arepo/MTNG-L500-4320-A/', type=str)
parser.add_argument('--snapnum', default=264, type=int)
args = parser.parse_args()

print('')
print(f'>>> Add profiles params to halo data <<<')
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
    NFWRs          = data['NFWRs']      
    rho_c          = data['rho_c'].value    # [Msun / (kpc)^3]
    print(NFWRs.shape)
    
    profile_params, fitted_densities, fitted_radius = [], [], []
    
    idx = 0
    for r, rho, R200, Rs in tqdm(zip(radial_bins, densities, halo_R_Mean200, NFWRs)):
        
        # Preprocess the data
        r, rho = filter_profile(r, rho) # Remove zero densities in the centre
        if args.sim.startswith('MTNG'):
            start_idx = np.where(r > 80)[0][0] # [kpc]
            r, rho = r[start_idx:], rho[start_idx:]
        
        sp = splashback_properties(tag=f'halo_{idx}', x=r, y=rho, mean=rho_c,
                           INNER='EINASTO', OUTER='MATTER',
                           Comparison=True, Norm=True, smooth=True)

        sp.fit_DK14()
        
        sp.compute_splashback_radius()
        sp.DK14_derivatives()

        sp.density_fit_plot()      

        params = sp.fit_params
        profile_params.append(params)
        
        # Substitute the 2nd parameter with the Rs
        sp.fit_params[1] = Rs
        fitted_rho = sp.DK14_function(sp.fit_params, INNER=sp.inner_pro, OUTER=sp.outer_pro, MODE='MODEL')
        fitted_densities.append(fitted_rho)
        fitted_radius.append(sp.x)
        print("          Rs:",Rs)
        idx += 1
 
    profile_params = np.array(profile_params)
    print(profile_params.shape)

    #  Append data
    data['profile_params']   = profile_params
    data['fitted_densities'] = np.array(fitted_densities)
    data['fitted_radius']    = np.array(fitted_radius)
    
    # Save the data
    print(data.keys())
    np.save(halos_dir+fname, data)
    print('Saved!')