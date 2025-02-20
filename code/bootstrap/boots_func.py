from func import *
import numpy as np
import math

def init_mass_bins(mass_1d_array):
    # Set up the mass bins
    log10_mass = np.log10(mass_1d_array) 
    min_log10_mass, max_log10_mass = np.min(log10_mass), np.max(log10_mass)
    # print(f'Min log10 mass: {min_log10_mass}, Max log10 mass: {max_log10_mass}')
    bin_start, bin_end, bin_width = math.floor(min_log10_mass), math.ceil(max_log10_mass), 1
    num_bins = int((bin_end-bin_start)/bin_width)
    mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun
    print('mass bins: ', mass_bins[:-1])
    return mass_bins
    
# Bootstrap setup
def bootstrap_all_features(args, params, data, formation_z=None):
    
    z, h, rho_c = params[0], params[1], params[2]
    radial_bins, densities, halo_M_Mean200, halo_R_Mean200 = data[0], data[1], data[2], data[3]
    total_num_halos = halo_M_Mean200.shape[0]
    
    mass_bins = init_mass_bins(halo_M_Mean200)
    
    # Initialize the results
    num_bins = len(mass_bins) - 1
    results = np.empty((num_bins, 6, args.Nboots))   # [med_mass, Rsp, depth, min_grad, width_dimless, width_physical]
    
    # Bootstrap
    valid_boots = 0
    while valid_boots < args.Nboots: 
        
        # Random selection of halos
        indices = np.random.randint(0, total_num_halos, args.Nsample)
        # Select densities and masses
        select_radii     = radial_bins[indices]    # [kpc]
        select_densities = densities[indices]      # [Msun / (kpc)^3]
        select_masses    = halo_M_Mean200[indices] # [10^10 Msun]
        select_r200      = halo_R_Mean200[indices] # [kpc]
        del indices
        
        # Calculating the number of halos in each cut
        num_halos_per_bin = count_halos(select_masses, 
                                        [10**m for m in mass_bins[:-1]], 
                                        [10**m for m in mass_bins[1:]])
        print('number of halos: ', num_halos_per_bin)
        if all(x > 2 for x in num_halos_per_bin):
            
            for i in range(num_bins):
                # Select halos in the cut and compute density profiles
                raw_profiles = stacked_density_profile(select_radii, select_densities, select_masses, select_r200, 
                                                    10**mass_bins[i], 10**mass_bins[i+1], rho_c)
                radius, rho, rho_err = raw_profiles[0][1:], raw_profiles[1][1:], raw_profiles[2][1:] # [dimensionless]
                _, R200_median = raw_profiles[3], raw_profiles[4]                             # [kpc]
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
                    print('This bootstrap is abandoned!')
                    results[:, :, valid_boots] = -1
                    break
                else:
                    # Fit the slope
                    fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho)) # [dimensionless]
                    
                    # Compute the median mass in the mass cut
                    med_mass = compute_median(select_masses, 10**mass_bins[i], 10**mass_bins[i+1])
                    
                    # Compute Rsp
                    physical_fitted_radius = fitted_radius * R200_median
                    Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
                    
                    # depth
                    min_grad = np.min(fitted_slope)
                    min_grad_idx = np.argmin(fitted_slope)
                    print(f'min grad index: {min_grad_idx}')
                    if min_grad_idx < 500:
                        print('This bootstrap is abandoned!')
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
                        results[i, :, valid_boots] = med_mass, Rsp, depth, min_grad, width_dimless, width
        
            ### Only the for loop is complete, update valid_boots
            # Updata counts
            if np.all(results[:, :, valid_boots] != -1):
                valid_boots += 1
                print(f'Nboots updated: {valid_boots}')
            else:
                print(f'Nboots not updated: {valid_boots} ')
            print('')
                
                
            
    # Get the statistical results 
    final_results = {'z': z, 'h': h, 'mass_bins': mass_bins[:-1]}
    if formation_z is not None:
        final_results['formation_z'] = formation_z
    final_results['med_mass']       = np.percentile(results[:, 0, :], [16, 50, 84], axis=1)
    final_results['Rsp']            = np.percentile(results[:, 1, :], [16, 50, 84], axis=1)
    final_results['depth']          = np.percentile(results[:, 2, :], [16, 50, 84], axis=1)
    final_results['abs_depth']      = np.percentile(results[:, 3, :], [16, 50, 84], axis=1)    
    final_results['width_dimless']  = np.percentile(results[:, 4, :], [16, 50, 84], axis=1)
    final_results['width_physical'] = np.percentile(results[:, 5, :], [16, 50, 84], axis=1)  
    final_results['full_results']   = results  
    return final_results