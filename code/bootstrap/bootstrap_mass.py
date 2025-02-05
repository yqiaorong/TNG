import os
import numpy as np
from func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1024, type=int)
# parser.add_argument('--bin_start',default=None, type=float)
# parser.add_argument('--bin_end',  default=None, type=float)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per mass cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Load halos data
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fname = os.listdir(halos_dir)[0]
print(halos_fname)
data = np.load(os.path.join(halos_dir, halos_fname), allow_pickle=True).item()

z            = data['z']
scale_factor = data['scale_factor']
h            = data['h']
rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 Msun]
densities      = data['densities']      # [Msun / (kpc)^3]
radial_bins    = data['radial_bins']    # [kpc]
del data

total_num_halos = halo_M_Mean200.shape[0]
print(f'total number of halos: {total_num_halos}')
print(halo_R_Mean200.shape, halo_M_Mean200.shape, densities.shape, radial_bins.shape)



# Make mass cuts
bin_start, bin_end, bin_width = int(halos_fname[4:6])/10, int(halos_fname[7:9])/10, 0.5
num_bins = int((bin_end-bin_start)/bin_width)
print(f'The number of mass bins: {num_bins}')
mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun
print(mass_bins[:-1])
print(f'The current mass range: 10^{bin_start+10} ~ 10^{bin_end+10} MSun')


    
# Bootstrap setup
Nsample, Nboots = args.Nsample, args.Nboots
results = np.empty((num_bins, 6, Nboots))   # [med_mass, Rsp, depth, min_grad, width_dimless, width_physical]

valid_boots = 0
while valid_boots < Nboots: 
    
    # Random selection of halos
    indices = np.random.randint(0, total_num_halos, Nsample)
    # Select densities and masses
    select_radii     = radial_bins[indices]    # [kpc]
    select_densities = densities[indices]      # [Msun / (kpc)^3]
    select_masses    = halo_M_Mean200[indices] # [10^10 Msun]
    select_r200      = halo_R_Mean200[indices] # [kpc]
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
            
            # Remove the radius < gravitational softening length ###
            # -----------------------------------------------------------------------------------------
            radius_phys = radius * R200_median # [kpc]
            # Get radius which is larger than the gravitational softening length
            if args.sim.startswith('MTNG'):
                rsoft = 80 # [kpc]
                start_idx = np.where(radius_phys > rsoft)[0][0]
                radius, rho, rho_err = radius[start_idx:], rho[start_idx:], rho_err[start_idx:]
                print(start_idx)
                
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
                
                # Compute the median mass in the mass cut
                med_mass = compute_median_mass(select_masses, mass_bins[i])
                
                # Compute Rsp
                physical_fitted_radius = fitted_radius * R200_median
                Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
                
                # depth
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
                results[i, :, valid_boots] = med_mass, Rsp, depth, min_grad, width_dimless, width
       
        ### Only the for loop is complete, update valid_boots
        # Updata counts
        valid_boots += 1
        print(f'Nboots: {valid_boots}')
        print('')
            
            
        
# Get the statistical results 
final_results = {'z': z, 'h': h, 'mass_bins': mass_bins[:-1]}
final_results['med_mass']       = np.percentile(results[:, 0, :], [16, 50, 84], axis=1)
final_results['Rsp']            = np.percentile(results[:, 1, :], [16, 50, 84], axis=1)
final_results['depth']          = np.percentile(results[:, 2, :], [16, 50, 84], axis=1)
final_results['abs_depth']      = np.percentile(results[:, 3, :], [16, 50, 84], axis=1)    
final_results['width_dimless']  = np.percentile(results[:, 4, :], [16, 50, 84], axis=1)
final_results['width_physical'] = np.percentile(results[:, 5, :], [16, 50, 84], axis=1)  
final_results['full_results']   = results    

# print(f'final_results shape (bin, type, percentile): {final_results.shape}')
# for i in range(num_bins):
#     print(f'bin {i}: ')
#     print(f'Rsp: {final_results[i, 0]}')
#     print(f'depth: {final_results[i, 1]}')
#     print(f'width dimless: {final_results[i, 2]}')
#     print(f'width physical: {final_results[i, 3]}')
#     print(f'depth absolute: {final_results[i, 4]}')
#     print('')
    
# # Check the index of median value
# origin_indices = []
# for i, cut in enumerate(results):
#     indices = np.argsort(cut[0,:]) # Rsp
#     origin_idx = np.where(indices == int(Nboots/2))[0][0]
#     print(f'cut {mass_bins[i]}: median boots idx = {origin_idx}')
#     origin_indices.append(origin_idx)
    
    
    
# Save the results
save_stats_dir = f'result/bootstrap_stats/with_mass/{args.sim}/Nboots_{Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)
    
# save_data = {'z': z, 'h': h, 'mass_bins': mass_bins[:-1],
#              'full_results': results, 
#              'final_results': final_results, 
#              'median_idx_in_boots': origin_indices}
np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', final_results)