import os
import numpy as np
from boots_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--profile',    default='DK14', type=str)
parser.add_argument('--sim',        default='MTNG/Hydro-Arepo/MTNG-L500-4320-A/', type=str)
parser.add_argument('--Nsample',    default=10000,type=int)
parser.add_argument('--Nboots',     default=1024, type=int)
parser.add_argument('--x_type',     default=None, type=str) 
parser.add_argument('--mass_start', default=None, type=float)
parser.add_argument('--mass_end',   default=None, type=float)
parser.add_argument('--reject_limit',default=1000, type=int)
args = parser.parse_args()
args.snapnum = 264

print('')
print(f'>>> Bootstrap splashback features per {args.x_type} cuts per mass bin <<<')
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
    halo_M_Mean200 = data['halo_M_Mean200'] # [Msun]
    print(halo_M_Mean200)
    x_data = data[args.x_type]  
    densities    = data['densities']      # [Msun / (kpc)^3]
    radial_bins  = data['radial_bins']    # [kpc] (num halos, num radial bins,)
    del data
    print(f'Original data shape: {x_data.shape}')
    
    # Select a subset of halos if bin_data contains NaNs
    valid_indices = ~np.isnan(x_data)    
    # If args.bin_type is 'NFWconc', constrain the value between 1 and 200
    if args.x_type == 'NFWconc':
        valid_indices &= (x_data >= 1) & (x_data <= 200) 
    elif args.x_type == 'accretions':
        valid_indices &= (x_data >= 0) & (x_data <= 6) 
        
    radial_bins    = radial_bins[valid_indices, :]  
    densities      = densities[valid_indices]
    x_data         = x_data[valid_indices]
    halo_R_Mean200 = halo_R_Mean200[valid_indices]
    halo_M_Mean200 = halo_M_Mean200[valid_indices]
    
    # Select data within a mass cut
    mass_idx = np.where((halo_M_Mean200 >= 10**args.mass_start) & (halo_M_Mean200 < 10**args.mass_end))[0]
    radial_bins    = radial_bins[mass_idx, :]
    densities      = densities[mass_idx]
    x_data         = x_data[mass_idx]
    halo_R_Mean200 = halo_R_Mean200[mass_idx]
    halo_M_Mean200 = halo_M_Mean200[mass_idx]
    print(f'Valid data shape: {x_data.shape}')



if args.x_type == 'accretions':
    x_start, x_end, x_width = 1, 6, 0.5
elif args.x_type == 'NFWconc':
    x_start, x_end, x_width = 0, 40, 2
elif args.x_type == 'peakHeight':
    x_start, x_end, x_width = 0, 6, 0.2
elif args.x_type in ['formz', 'formzOLD', 'formzSub']:
    x_start, x_end, x_width = 0, 1.6, 0.2
elif args.x_type == 'mergerz':
    x_start, x_end, x_width = 0, 10, 1
else:
    x_start, x_end, x_width = None, None, None
    
    

if x_data.shape[0] == 0:
    print('No valid data found!')
    exit()
else:
    # Bootstrap setup
    final_results = bootstrap_all_features(args, 
                                            [z, h, rho_c], 
                                            [radial_bins, densities, x_data, halo_R_Mean200],
                                            args.x_type, plot_x_type_name=args.x_type+f'_mass_{str(int(args.mass_start*10))}',                                     
                                            x_start=x_start, x_end=x_end, x_width=x_width                             
                                            ) 
        
    # Save the results
    save_stats_dir = f'result/bootstrap_stats_{args.profile}/with_{args.x_type}_perMassCut/{args.sim}/Nboots_{args.Nboots}/'
    if not os.path.exists(save_stats_dir):
        os.makedirs(save_stats_dir)

    np.save(save_stats_dir+f'snap_{args.snapnum}_mass_{str(int(args.mass_start*10))}_Rsp_stats', final_results)
    print('Saved!')