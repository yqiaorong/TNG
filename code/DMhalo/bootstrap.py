import h5py
import os
import numpy as np
from func import *
import illustris_python as il
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--snapnum',default=99,     type=int)
parser.add_argument('--Nsample',default=10000,  type=int)
parser.add_argument('--bin_end',default=None,   type=float)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap Rsp <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Make mass cuts

# if args.snapnum == 8:
#     bin_end = 2
# elif args.snapnum == 13:
#     bin_end = 2.5
# elif args.snapnum == 17 or 21 or 25:
#     bin_end  = 3
# else:
#     bin_end = 1, 4

bin_start, bin_width = 1, 0.5
bin_end = args.bin_end
num_bins = int((bin_end-bin_start)/bin_width)
mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun/h
print(f'The current mass range: 10^{bin_start+10} ~ 10^{bin_end+10} MSun/h')

    
snap = args.snapnum

# Load redshift values
with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    z = 1 / scale_factor - 1
    h = header['HubbleParam']
    
    

### Load halos ###

# Halos data root dir
halos_dir = f'result/{root_dir}/sim_{boxsize}_{res}/snap_{snap}/densities'

# First round of rough selection of halos based on M200
Group_M_Mean200 = il.groupcat.loadHalos(basePath, snap, fields='Group_M_Mean200')
subset_idx = np.where((Group_M_Mean200 >= 10**bin_start) & (Group_M_Mean200 < 10**bin_end))[0]
halos_list = [f'halo_{idx}.npy' for idx in subset_idx]
print(f'The total halo numbers: {len(halos_list)}')
del Group_M_Mean200

# Save the result
save_data_dir = f'result/bootstrap/sim_{boxsize}_{res}'
if not os.path.exists(save_data_dir):
    os.makedirs(save_data_dir)
save_stats_dir = f'result/bootstrap_stats/sim_{boxsize}_{res}'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)

### Bootstrap ###

# Bootstrap setup
Nsample, Nboots = args.Nsample, 32
results = np.empty((num_bins, 3, Nboots))

valid_boots = 0
while valid_boots < Nboots: 
    
    # Select halos
    indices = np.random.randint(0, len(halos_list), Nsample)
    boots_halo_list = np.array(halos_list)[indices]
    print(f'Cross check: the number of selected halos in bootstrap: {len(boots_halo_list)}')
    del indices
    
    # Calculating the number of halos in each cut
    num_halos = count_halos(halos_dir, boots_halo_list, bin_start, bin_end)
    print(num_halos)
    
    if all(x > 1 for x in num_halos):
        
        for i in range(num_bins):
            # Select halos in the cut and compute density profiles
            raw_profiles = stacked_density_profile(halos_dir, boots_halo_list, [mass_bins[i], mass_bins[i+1]])
            radius, rho, rho_err = raw_profiles[0][1:], raw_profiles[1][1:], raw_profiles[2][1:] # [dimensionless]
            num_halo, R200_median = raw_profiles[3], raw_profiles[4]                             # [ckpc/h]
            del raw_profiles
            
            # Compute the slope
            radius, rho, rho_err = filter_profile(radius, rho, rho_err) # Remove zero densities in the centre
            slope = num_deriv(np.log(radius), np.log(rho))                                      # [dimensionless]
            slope_err = num_deriv_err(radius, rho, rho_err)                                     # [dimensionless]

            # Fit the density profiles
            fit_profiles = fit_profile_parametric(radius, rho, np.mean(rho_err, axis=1), 1)
            fitted_radius, fitted_rho = fit_profiles[0], fit_profiles[1]
            del fit_profiles
            
            # Fit the slope
            fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho))      # [dimensionless]

            # Plot the profile
            plot_profile(R200_median, radius, rho, rho_err, slope, slope_err, 
                        fitted_radius, fitted_rho, fitted_slope, 
                        [mass_bins[i], mass_bins[i+1]], num_halo, snap, 
                        f'{save_data_dir}/snap_{snap}', f'boots_{valid_boots}',
                        save_data=True)
            
            # Compute Rsp
            physical_fitted_radius = fitted_radius * R200_median * scale_factor / h
            Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
            
            # Rsp depth
            min_grad = np.min(fitted_slope)
            min_grad_idx = np.argmin(fitted_slope)
            print(f'min grad index: {min_grad_idx}')
            left_data = fitted_slope[:min_grad_idx]
            right_data = fitted_slope[min_grad_idx:]
            
            max_grad = np.max(right_data)
            depth = max_grad - min_grad
            
            # Width
            half_grad = min_grad + depth/2
            
            left_idx = np.argmin(np.abs(left_data - half_grad))
            right_idx = min_grad_idx + np.argmin(np.abs(right_data - half_grad))
            
            width = physical_fitted_radius[right_idx] - physical_fitted_radius[left_idx]
            
            # Append results
            results[i, :, valid_boots] = Rsp, depth, width
        del boots_halo_list
        
        # Updata counts
        valid_boots += 1
        print(f'Nboots: {valid_boots}')
        print('')
            
final_results = np.percentile(results, [16, 50, 84], axis=2).transpose(1,2,0)
print(f'final_results shape (bin, type, percentile): {final_results.shape}')
for i in range(num_bins):
    print(f'bin {i}: ')
    print(f'Rsp: {final_results[i, 0]}')
    print(f'depth: {final_results[i, 1]}')
    print(f'width: {final_results[i, 2]}')
    print('')
    
# Check the index of median value
origin_indices = []
for i, cut in enumerate(results):
    indices = np.argsort(cut[0,:]) # Rsp
    origin_idx = np.where(indices == int(Nboots/2))[0][0]
    print(f'cut {mass_bins[i]}: median boots idx = {origin_idx}')
    origin_indices.append(origin_idx)
    
save_data = {'z': z, 'h': h, 'start_mass_cut': bin_start, 
             'full_results': results, 
             'final_results': final_results, 
             'median_idx_in_boots': origin_indices}
np.save(save_stats_dir+f'/snap_{snap}_Rsp_stats', save_data)