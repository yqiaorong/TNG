import os
import numpy as np
from boots_func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default='MTNG/Hydro-Arepo/MTNG-L500-4320-A/', type=str)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1024, type=int)
parser.add_argument('--x_type',   default=None, type=str) # [ conc / peakHeight / accretions ]
parser.add_argument('--bin_type', default=None, type=str) # [ conc / peakHeight / accretions ]
parser.add_argument('--bin_start', default=None, type=float)
parser.add_argument('--bin_end',   default=None, type=float)
parser.add_argument('--reject_limit',default=2000, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features vs {args.x_type} per {args.bin_type} cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Load halos data
snaps = [129, 151, 179, 214, 237, 264]
all_z, all_h = [], []
all_R, all_x_data, all_bin_data, all_radii, all_densities = [], [], [], [], []
for snap in snaps:
    halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{snap}/final_densities/'
    halos_fnames = os.listdir(halos_dir)
    print('snap:', snap, halos_fnames)
    for fname in halos_fnames:
        data = np.load(halos_dir+fname, allow_pickle=True).item()
        
        all_z.append(data['z'])
        all_h.append(data['h'])

        halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
        x_data = data[args.x_type]
        if args.bin_type == 'mass':
            bin_data = data['halo_M_Mean200']   # [Msun]
        else:
            bin_data = data[args.bin_type]  
        densities    = data['densities'] / data['rho_c']  # [dimless]
        radial_bins  = data['radial_bins']    # [kpc] (num halos, num radial bins,)
        del data
        
        # Append data to lists
        all_R.append(halo_R_Mean200)
        all_x_data.append(x_data)
        all_bin_data.append(bin_data)
        all_radii.append(radial_bins)
        all_densities.append(densities)
del halo_R_Mean200, x_data, bin_data, densities, radial_bins

# Concatenate all dataW
all_R         = np.concatenate(all_R)
all_x_data    = np.concatenate(all_x_data)
all_bin_data  = np.concatenate(all_bin_data)
all_radii     = np.concatenate(all_radii)
all_densities = np.concatenate(all_densities)
print(all_R.shape, all_x_data.shape, all_bin_data.shape, all_radii.shape, all_densities.shape)  



# Select a subset of halos if bin_data contains NaNs
valid_indices = ~np.isnan(all_bin_data)    
# If args.bin_type is 'NFWconc', constrain the value between 1 and 200
if args.x_type == 'NFWconc':
    valid_indices &= (all_bin_data >= 1) & (all_bin_data <= 200) 
elif args.x_type == 'accretions':
    valid_indices &= (all_bin_data >= 0) & (all_bin_data <= 6) 
all_radii     = all_radii[valid_indices]  
all_densities = all_densities[valid_indices]
all_bin_data  = all_bin_data[valid_indices]
all_x_data    = all_x_data[valid_indices]
all_R         = all_R[valid_indices]



# Select data within a cut
if args.bin_type == 'mass':
    select_idx = np.where((all_bin_data >= 10**args.bin_start) & (all_bin_data < 10**args.bin_end))[0] 
else:
    select_idx = np.where((all_bin_data >= args.bin_start) & (all_bin_data < 10**args.bin_end))[0]
all_radii     = all_radii[select_idx]  
all_densities = all_densities[select_idx]
all_bin_data  = all_bin_data[select_idx]
all_x_data    = all_x_data[select_idx]
all_R         = all_R[select_idx]



# Shuffle the data based on all_x_data
shuffled_indices = np.random.permutation(len(all_x_data))
all_R         = all_R[shuffled_indices]
all_x_data    = all_x_data[shuffled_indices]
all_bin_data  = all_bin_data[shuffled_indices]
all_radii     = all_radii[shuffled_indices]
all_densities = all_densities[shuffled_indices]
print(all_R.shape, all_x_data.shape, all_bin_data.shape, all_radii.shape, all_densities.shape)  



# if args.bin_type == 'accretions':
#     bin_start, bin_end, bin_width = 1, 6, 0.5
# elif args.bin_type == 'NFWconc':
#     bin_start, bin_end, bin_width = 0, 40, 2
if args.x_type == 'peakHeight':
    x_start, x_end, x_width = 0, 6, 0.2
# elif args.bin_type in ['formz', 'formzOLD', 'formzSub']:
#     bin_start, bin_end, bin_width = 0, 1.6, 0.2
# elif args.bin_type == 'mergerz':
#     bin_start, bin_end, bin_width = 0, 10, 1
else:
    x_start, x_end, x_width = None, None, None

    

if all_bin_data.shape[0] == 0:
    print('No valid data found!')
    exit()
else:
    # Bootstrap setup
    final_results, _ = bootstrap_all_features(args, 
                                            [all_z, all_h], 
                                            [all_radii, all_densities, all_x_data, all_R],
                                            args.x_type, plot_x_type_name=args.x_type+f'_{args.bin_type}_{str(int(args.bin_start))}',                                      
                                            x_start=x_start, x_end=x_end, x_width=x_width                             
                                            ) 
        
    # Save the results
    save_stats_dir = f'result/bootstrap_stats/with_{args.x_type}_per{args.bin_type}Cut/{args.sim}/Nboots_{args.Nboots}/'
    if not os.path.exists(save_stats_dir):
        os.makedirs(save_stats_dir)

    np.save(save_stats_dir+f'{args.bin_type}_{str(args.bin_start)}_Rsp_stats', final_results)
    print('Saved!')