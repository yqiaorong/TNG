import os
import numpy as np
from boots_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',        default='MTNG/Hydro-Arepo/MTNG-L500-4320-A/', type=str)
parser.add_argument('--snapnum',    default=264,  type=int)
parser.add_argument('--Nsample',    default=5000, type=int)
parser.add_argument('--Nboots',     default=1024, type=int)
parser.add_argument('--bin_type',   default=None, type=str) 
parser.add_argument('--mass_start', default=None, type=float)
parser.add_argument('--mass_end',   default=None, type=float)
parser.add_argument('--reject_limit',default=1000, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per {args.bin_type} cuts per mass bin <<<')
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
    halo_M_Mean200 = halo_M_Mean200[valid_indices]
    
    # Select data within a mass cut
    mass_idx = np.where((halo_M_Mean200 >= 10**args.mass_start) & (halo_M_Mean200 < 10**args.mass_end))[0]
    radial_bins    = radial_bins[mass_idx, :]
    densities      = densities[mass_idx]
    bin_data       = bin_data[mass_idx]
    halo_R_Mean200 = halo_R_Mean200[mass_idx]
    halo_M_Mean200 = halo_M_Mean200[mass_idx]
    print(f'Valid data shape: {bin_data.shape}')



# NGWconc: 5; peakHeight: 0.2; accretions: 0.5.
if args.bin_type == 'accretions':
    bin_start, bin_end, bin_width = 1, 6, 0.5
elif args.bin_type == 'NFWconc':
    bin_start, bin_end, bin_width = 0, 40, 2
elif args.bin_type == 'peakHeight':
    bin_start, bin_end, bin_width = 0, 6, 0.2
elif args.bin_type in ['formz', 'formzOLD', 'formzSub']:
    bin_start, bin_end, bin_width = 0, 1.6, 0.2
elif args.bin_type == 'mergerz':
    bin_start, bin_end, bin_width = 0, 10, 1
else:
    bin_start, bin_end, bin_width = None, None, None
    
    

if bin_data.shape[0] == 0:
    print('No valid data found!')
    exit()
else:
    # Bootstrap setup
    final_results, _ = bootstrap_all_features(args, 
                                            [z, h, rho_c], 
                                            [radial_bins, densities, bin_data, halo_R_Mean200],
                                            args.bin_type, plot_bin_type_name=args.bin_type+f'_mass_{str(int(args.mass_start*10))}',                                     
                                            bin_start=bin_start, bin_end=bin_end, bin_width=bin_width                             
                                            ) 
        
    # Save the results
    save_stats_dir = f'result/bootstrap_stats/with_{args.bin_type}_perMassCut/{args.sim}/Nboots_{args.Nboots}/'
    if not os.path.exists(save_stats_dir):
        os.makedirs(save_stats_dir)

    np.save(save_stats_dir+f'snap_{args.snapnum}_mass_{str(int(args.mass_start*10))}_Rsp_stats', final_results)
    print('Saved!')
    
    # # Save the parameters
    # save_stats_dir2 = f'result/bootstrap_stats2/with_{args.bin_type}_perMassCut/{args.sim}/Nboots_{args.Nboots}/'
    # if not os.path.exists(save_stats_dir2):
    #     os.makedirs(save_stats_dir2)

    # np.save(save_stats_dir2+f'snap_{args.snapnum}_mass_{str(int(args.mass_start*10))}_Rsp_stats', final_params_results)
    # print('Saved!')