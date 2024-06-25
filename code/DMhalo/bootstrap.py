import h5py
import os
import numpy as np
from matplotlib import pyplot as plt
from func import *
import illustris_python as il

root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Make mass cuts
bin_start, bin_end, bin_width = 2.5, 4, 0.5
num_bins = int((bin_end-bin_start)/bin_width)
mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun/h


    
snap = 33
print(f'The current snapshot: {snap}')

# Load redshift values
with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
    header = dict(f['Header'].attrs.items())
    scale_factor = header['Time']
    z = 1 / scale_factor - 1
    h = header['HubbleParam']
    
    

# Load halos 
halos_dir = f'result/{root_dir}/sim_{boxsize}_{res}/snap_{snap}/densities'
halos_list = [os.path.join(halos_dir, fname) for fname in os.listdir(halos_dir)]
print(f'The total halo numbers: {len(halos_list)}')

# First round of rough selection of halos
halos_list = select_halos(halos_list, bin_start, bin_end)
print(f'The total halo numbers: {len(halos_list)}')



# Bootstrap setup
Nsample, Nboots = 1000, 32
results = np.empty((num_bins, 3, Nboots))

# Bootstrap
valid_boots = 0
while valid_boots < Nboots: 
    
    # Select halos
    indices = np.random.randint(0, len(halos_list), Nsample)
    sub_list = np.array(halos_list)[indices]
    print(f'Cross check: the number of selected halos: {len(sub_list)}')
    
    # Calculating the number of halos in each cut
    num_halos = count_halos(sub_list, bin_start, bin_end)
    print(num_halos)
    
    if all(x > 10 for x in num_halos):
        
        for i in range(num_bins):
            # Select halos in the cut and compute density profiles
            raw_profiles = stacked_density_profile(sub_list, [mass_bins[i], mass_bins[i+1]])
            radius, rho, rho_err = raw_profiles[0][1:], raw_profiles[1][1:], raw_profiles[2][1:] # [dimensionless]
            num_halo, R200_median = raw_profiles[3], raw_profiles[4]                             # [ckpc/h]
            del raw_profiles
            
            # Compute the slope
            slope = num_deriv(np.log(radius), np.log(rho))                            # [dimensionless]
            slope_err = num_deriv_err(radius, rho, rho_err)
            
            # Fit the density profiles
            fit_profiles = fit_profile_parametric(radius, rho, np.mean(rho_err, axis=1), 1)
            fitted_radius, fitted_rho = fit_profiles[0], fit_profiles[1]
            del fit_profiles
            
            # Fit the slope
            fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho))      # [dimensionless]

            # Plot the profile
            plot_profile(radius, rho, rho_err, slope, slope_err, 
                        fitted_radius, fitted_rho, fitted_slope, 
                        [mass_bins[i], mass_bins[i+1]], num_halo, snap, 
                        f'result/bootstrap/snap_{snap}')
            
            
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
            
            width = fitted_radius[right_idx] - fitted_radius[left_idx]
            
            # Append results
            results[i, :, valid_boots] = Rsp, depth, width
            
        # Updata counts
        valid_boots += 1
        print(f'Nboots: {valid_boots}')
        print('')
            
final_results = np.percentile(results, [16, 50, 84], axis=2)
print(f'final_results shape (percentile, bin, type)')
for i in range(int(len(mass_bins)-1)):
    print(f'bin {i}: ')
    print(f'Rsp: {final_results[:, i, 0]}')
    print(f'depth: {final_results[:, i, 1]}')
    print(f'depth: {final_results[:, i, 2]}')
    print('')
    
# Save the result
save_dir = f'result/bootstrap/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

np.save(save_dir+f'snap_{snap}_Rsp_stats', final_results)