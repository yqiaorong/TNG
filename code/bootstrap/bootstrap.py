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
parser.add_argument('--bin_type', default=None, type=str) # [ conc / peakHeight / accretions ]
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per {args.bin_type} cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load halos data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    print(data.keys())
    z            = data['z']
    scale_factor = data['scale_factor']
    h            = data['h']
    rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
    halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
    bin_data = data[args.bin_type]  
    densities    = data['densities']      # [Msun / (kpc)^3]
    radial_bins  = data['radial_bins']    # [kpc] (num halos, num radial bins,)
    del data
    print(f'Original data shape: {bin_data.shape}')
    
    # Select a subset of halos if bin_data contains NaNs
    valid_indices = ~np.isnan(bin_data)    
    # If args.bin_type is 'NFWconc', constrain the value between 1 and 200
    if args.bin_type == 'NFWconc':
        valid_indices &= (bin_data >= 1) & (bin_data <= 200) 
    elif args.bin_type in ['accretions', 'accretionsOLD']:
        valid_indices &= (bin_data >= 0) & (bin_data <= 6) 
        
    radial_bins    = radial_bins[valid_indices, :]  
    densities      = densities[valid_indices]
    bin_data       = bin_data[valid_indices]
    halo_R_Mean200 = halo_R_Mean200[valid_indices]
    print(f'Valid data shape: {bin_data.shape}')
    
if args.bin_type in ['mergerz', 'formz']:
    bins = np.unique(bin_data)
    if args.bin_type == 'formz':
        bins = bins[:-1]
else:
    bins = None

# Bootstrap setup
final_results = bootstrap_all_features(args, 
                                [z, h, rho_c], 
                                [radial_bins, densities, bin_data, halo_R_Mean200],
                                args.bin_type, 
                                # bins=bins,                                        # For mergerz, formz only
                                bin_width = 0.5                                   # NGWconc: 10; peakHeight: 0.2; accretions: 0.5.
                                ) 
    
# Save the results
save_stats_dir = f'result/bootstrap_stats/with_{args.bin_type}/{args.sim}/Nboots_{args.Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)

np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', final_results)