import numpy as np
import os
import argparse
from boots_func import *


# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per mass cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load the data
# ---------------------------------------------------------------------
load_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
fname = os.listdir(load_dir)[0]
data = np.load(load_dir + fname, allow_pickle=True).item()

halo_M      = data['halo_M_Mean200'] 
halo_R      = data['halo_R_Mean200']
halo_radii  = data['radial_bins']
halo_densities = data['densities']
formation_z = data['formation_time']

round_z  = np.round(formation_z, 1)
unique_z = np.unique(round_z)
# remove nan from unique_z
unique_z = unique_z[~np.isnan(unique_z)]
del formation_z


# Bin the halos both by formation redshift
# ---------------------------------------------------------------------
for form_z in unique_z:
    print(f'formation z: {form_z}')
    
    # Select the subset of data
    z_idx = np.where(round_z == form_z)[0]
    select_M200 = halo_M[z_idx]
    select_R200 = halo_R[z_idx]
    select_radii = halo_radii[z_idx]
    select_densities = halo_densities[z_idx]
    print(f'Number of halos: {z_idx.shape[0]}')
    
    if z_idx.shape[0] < 10:
        print(f'{z_idx.shape[0]} halos: Not enough for this formation redshift')
    else:
        # Bootstrap by mass cut
        results_at_form_z = bootstrap_all_features(args, 
                                                [data['z'], data['h'], data['rho_c']], 
                                                [select_radii, select_densities, select_M200, select_R200], 
                                                formation_z=form_z)

        # Save the results
        if results_at_form_z != {}:
            save_stats_dir = f'result/bootstrap_stats/with_formation_z/{args.sim}/Nboots_{args.Nboots}/snap_{args.snapnum}/'
            if not os.path.exists(save_stats_dir):
                os.makedirs(save_stats_dir)

            np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats_at_{form_z}', results_at_form_z)
    print('')