import os
import math
import numpy as np
from func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',      default=None, type=str)
parser.add_argument('--snapnum',  default=None, type=int)
parser.add_argument('--Nsample',  default=10000,type=int)
parser.add_argument('--Nboots',   default=1024, type=int)
parser.add_argument('--accret_start',default=None, type=int)
parser.add_argument('--accret_end',  default=None, type=int)
args = parser.parse_args()

print('')
print(f'>>> Bootstrap splashback features per accretion rate cuts <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')



# Load halos data
# -----------------------------------------------------------------------------------------
halos_dir = f'result/DMhalo_density_profiles/{args.sim}/snap_{args.snapnum}/final_densities/'
halos_fname = os.listdir(halos_dir)[0]
print(halos_fname)
data = np.load(os.path.join(halos_dir, halos_fname), allow_pickle=True).item()

z            = data['z']
scale_factor = data['scale_factor']
h            = data['h']
rho_c        = data['rho_c']            # [(Msun) / (kpc)^3]
halo_R_Mean200 = data['halo_R_Mean200'] # [kpc]
accretion_rate = data['accretion_rate'] # [dimless]
densities      = data['densities']      # [Msun / (kpc)^3]
radial_bins    = data['radial_bins']    # [kpc]
del data

# Check the accretion rates
print(accretion_rate.shape, halo_R_Mean200.shape, densities.shape, radial_bins.shape)

# Filter out nan values
valid_indices  = np.where(~np.isnan(accretion_rate))[0]
halo_R_Mean200 = halo_R_Mean200[valid_indices]
accretion_rate = accretion_rate[valid_indices]
densities      = densities[valid_indices]
radial_bins    = radial_bins[valid_indices]
print(accretion_rate.shape, halo_R_Mean200.shape, densities.shape, radial_bins.shape)

total_num_halos = accretion_rate.shape[0]



# Create accretion rate bins
# -----------------------------------------------------------------------------------------
if accretion_rate.shape[0] == 0:
	exit()
else:
	print(min(accretion_rate), max(accretion_rate))
	# Bin the accretion rate cuts by 1
	# bin_start, bin_end, bin_width = math.floor(min(accretion_rate)), math.ceil(max(accretion_rate)), 1
	bin_start, bin_end, bin_width = args.accret_start, args.accret_end, 1
	num_bins = int((bin_end-bin_start)/bin_width)
	print(f'The number of acretion rate bins: {num_bins}')
	accret_bins = np.arange(bin_start, bin_end+bin_width, bin_width)
	print(accret_bins)
	print(f'The current accretion rate range: {bin_start} ~ {bin_end}')



# Bootstrap setup
# -----------------------------------------------------------------------------------------
Nsample, Nboots = args.Nsample, args.Nboots
results = np.empty((num_bins, 6, Nboots))   # [med_accret, Rsp, depth, min_grad, width_dimless, width_physical]

valid_boots = 0
while valid_boots < Nboots: 
    
    print('Current boots: ', valid_boots)
    # Random selection of halos
    indices = np.random.randint(0, total_num_halos, Nsample)
    # Select densities and masses
    select_radii     = radial_bins[indices]    # [kpc]
    select_densities = densities[indices]      # [Msun / (kpc)^3]
    select_accrets   = accretion_rate[indices] 
    select_r200      = halo_R_Mean200[indices] # [kpc]
    del indices
    
    # Calculating the number of halos in each cut
    num_halos_per_bin = count_halos(select_accrets, accret_bins[:-1], accret_bins[1:])
    print(num_halos_per_bin)
    
    if all(x > 2 for x in num_halos_per_bin):
        
        for i in range(num_bins):
            # Select halos in the cut and compute density profiles
            raw_profiles = stacked_density_profile(select_radii, select_densities, select_accrets, select_r200, 
                                                   accret_bins[i], accret_bins[i+1], rho_c)
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
                results[:, :, valid_boots] = -1
                break
            else:
                # Fit the slope
                fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho)) # [dimensionless]
                
                # Compute the median mass in the mass cut
                med_accret = compute_median(select_accrets, accret_bins[i], accret_bins[i+1])
                
                # Compute Rsp
                physical_fitted_radius = fitted_radius * R200_median
                Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
                
                # depth
                min_grad = np.min(fitted_slope)
                min_grad_idx = np.argmin(fitted_slope)
                print(f'min grad index: {min_grad_idx}')
                if min_grad_idx < 100:
                    print('This bootstrap is abandoned! ')
                    results[:, :, valid_boots] = -1
                    break
                else:
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
                    results[i, :, valid_boots] = med_accret, Rsp, depth, min_grad, width_dimless, width
                    print(med_accret, Rsp, depth, min_grad, width_dimless, width)
                    
        ### Only the for loop is complete, update valid_boots
        # Updata counts
        if np.all(results[:, :, valid_boots] != -1):
            valid_boots += 1
            print(f'Nboots updated: {valid_boots}')
        else:
            print(f'Nboots not updated: {valid_boots} ')
        print('')
            
            
        
# Get the statistical results 
final_results = {'z': z, 'h': h, 'accret_bins': accret_bins[:-1]}
final_results['med_accret']     = np.percentile(results[:, 0, :], [16, 50, 84], axis=1)
final_results['Rsp']            = np.percentile(results[:, 1, :], [16, 50, 84], axis=1)
final_results['depth']          = np.percentile(results[:, 2, :], [16, 50, 84], axis=1)
final_results['abs_depth']      = np.percentile(results[:, 3, :], [16, 50, 84], axis=1)    
final_results['width_dimless']  = np.percentile(results[:, 4, :], [16, 50, 84], axis=1)
final_results['width_physical'] = np.percentile(results[:, 5, :], [16, 50, 84], axis=1)  
final_results['full_results']   = results    
    
    
    
# Save the results
save_stats_dir = f'result/bootstrap_stats/with_accret/{args.sim}/Nboots_{Nboots}/'
if not os.path.exists(save_stats_dir):
    os.makedirs(save_stats_dir)
    
np.save(save_stats_dir+f'snap_{args.snapnum}_Rsp_stats', final_results)