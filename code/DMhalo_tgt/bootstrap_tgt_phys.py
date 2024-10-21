import h5py
import os
import numpy as np
from func import *
import illustris_python as il
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1024, type=int)
parser.add_argument('--bin_start',default=None, type=float)
parser.add_argument('--bin_end',  default=None, type=float)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap Physical Rsp <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Load redshift values (Alternative!!!)
halos_dir = f'result/DMhalo_density_profiles_phys/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_list = os.listdir(halos_dir)
sample_file = halos_list[0]
sample_data = np.load(os.path.join(halos_dir, sample_file), allow_pickle=True).item()
print(sample_data.keys())
z = sample_data['z']
scale_factor = sample_data['scale_factor']
h = sample_data['h']
rho_c = sample_data['rho_c'] # [(Msun) / (kpc)^3]
del sample_data



# Make mass cuts
bin_start, bin_end, bin_width = args.bin_start, args.bin_end, 0.5
num_bins = int((bin_end-bin_start)/bin_width)
print(f'The number of mass bins: {num_bins}')
mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun
print(mass_bins[:-1])
print(f'The current mass range: 10^{bin_start+10} ~ 10^{bin_end+10} MSun')

# Halos data root dir

print(halos_list)

# Load halos data
total_num_halos = 0
halo_R_Mean200, halo_M_Mean200, densities, radial_bins = [], [], [], []
for ifname, fname in enumerate(halos_list):
    data = np.load(f'{halos_dir}/{fname}', allow_pickle=True).item()

    if ifname == 0:
        halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
        halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 Msun]
        densities = data['densities']           # [Msun / (kpc)^3]
        radial_bins = data['radial_bins']       # [kpc]
    else:
        halo_R_Mean200 = np.concatenate((halo_R_Mean200, data['halo_R_Mean200']), axis=0)
        halo_M_Mean200 = np.concatenate((halo_M_Mean200, data['halo_M_Mean200']), axis=0)
        densities = np.concatenate((densities, data['densities']), axis=0)
        radial_bins = np.concatenate((radial_bins, data['radial_bins']), axis=0)
    
    total_num_halos = total_num_halos + data['halo_R_Mean200'].shape[0]
    del data
    
print(halo_R_Mean200.shape, halo_M_Mean200.shape, densities.shape, radial_bins.shape)
print(f'total number of halos: {total_num_halos}')

# Save plot dir 
save_data_dir = f'result/bootstrap_phys/{args.sim}/snap_{args.snapnum}/Nboots_{args.Nboots}/'
if not os.path.exists(save_data_dir):
    os.makedirs(save_data_dir)
    
# Bootstrap setup
Nsample, Nboots = args.Nsample, args.Nboots
results = np.empty((num_bins, 4, Nboots))

valid_boots = 0
while valid_boots < Nboots: 
    
    # Random selection of halos
    indices = np.random.randint(0, total_num_halos, Nsample)
    # Select densities and masses
    select_radii = radial_bins[indices]     # [kpc]
    select_densities = densities[indices]   # [Msun / (kpc)^3]
    select_masses = halo_M_Mean200[indices] # [10^10 Msun]
    select_r200 = halo_R_Mean200[indices]   # [kpc]

    del indices
    
    # Calculating the number of halos in each cut
    num_halos_per_bin = count_halos_based_on_mass(select_masses, mass_bins[:-1])
    print(num_halos_per_bin)
    if all(x > 1 for x in num_halos_per_bin):
        
        for i in range(num_bins):
            # Select halos in the cut and compute density profiles
            raw_profiles = stacked_density_profile(select_radii, select_densities, select_masses, select_r200, 
                                                   mass_bins[i], h, scale_factor, rho_c)
            radius, rho, rho_err = raw_profiles[0][1:], raw_profiles[1][1:], raw_profiles[2][1:] # [dimensionless]
            num_halo, R200_median = raw_profiles[3], raw_profiles[4]                             # [kpc]
            del raw_profiles

            # Compute the slope
            radius, rho, rho_err = filter_profile(radius, rho, rho_err) # Remove zero densities in the centre
            slope = num_deriv(np.log(radius), np.log(rho))                                      # [dimensionless]
            slope_err = num_deriv_err(radius, rho, rho_err)                                     # [dimensionless]
            
            # Fit the density profiles
            fit_profiles = fit_profile_parametric(radius, rho, np.mean(rho_err, axis=1), 1)
            fitted_radius, fitted_rho = fit_profiles[0], fit_profiles[1]
            del fit_profiles
            
            ### If the optimal params are not found! ###
            if np.all(fitted_rho) == 0:
                print('This bootstrap is abandoned! ')
                break
            else:
                # Fit the slope
                fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho)) # [dimensionless]

                # Plot the profile
                plot_profile(R200_median, radius, rho, rho_err, slope, slope_err, 
                            fitted_radius, fitted_rho, fitted_slope, 
                            [mass_bins[i], mass_bins[i+1]], num_halo, args.snapnum, 
                            save_data_dir, f'boots_{valid_boots}',
                            save_data=True)
                
                # Compute Rsp
                physical_fitted_radius = fitted_radius * R200_median
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
                width_dimless = fitted_radius[right_idx] - fitted_radius[left_idx]
                
                # Append results
                results[i, :, valid_boots] = Rsp, depth, width_dimless, width
       
        ### Only the for loop is complete, update valid_boots
        else:
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
    print(f'width dimless: {final_results[i, 2]}')
    print(f'width physical: {final_results[i, 3]}')
    print('')
    
# Check the index of median value
origin_indices = []
for i, cut in enumerate(results):
    indices = np.argsort(cut[0,:]) # Rsp
    origin_idx = np.where(indices == int(Nboots/2))[0][0]
    print(f'cut {mass_bins[i]}: median boots idx = {origin_idx}')
    origin_indices.append(origin_idx)
    
    
    
# Save the result
save_stats_dir = f'result/bootstrap_stats_phys/{args.sim}/Nboots_{Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)
    
save_data = {'z': z, 'h': h, 'mass_bins': mass_bins[:-1],
             'full_results': results, 
             'final_results': final_results, 
             'median_idx_in_boots': origin_indices}
np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', save_data)