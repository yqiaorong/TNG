from func import *
import numpy as np
import math

def init_mass_bins(mass_1d_array, bin_width):
    
    # Set up the mass bins
    log10_mass = np.log10(mass_1d_array) 
    min_log10_mass, max_log10_mass = np.min(log10_mass), np.max(log10_mass)
    # print(f'Min log10 mass: {min_log10_mass}, Max log10 mass: {max_log10_mass}')
    bin_start, bin_end = math.floor(min_log10_mass*2)/2, math.ceil(max_log10_mass*2)/2
    mass_bins = np.arange(bin_start, bin_end+bin_width, bin_width)
    print('mass bins: ', mass_bins)
    
    # Check the number of halos in the largest bin. if it's less than 2, reject this bin
    new_mass_bins = []
    for mass_bin in mass_bins[:-1]:
        num_halos_in_bin = np.where((log10_mass >= mass_bin) & (log10_mass < mass_bin+bin_width))[0].shape[0]
        # print(mass_bin, num_halos_in_bin)
        if num_halos_in_bin > 100:
            new_mass_bins.append(mass_bin)
            
    if len(new_mass_bins) != 0:
        new_mass_bins.append(new_mass_bins[-1]+bin_width)
    del mass_bins
        
    print('mass bins: ', new_mass_bins)

    return new_mass_bins
    
# Bootstrap setup
def bootstrap_all_features(args, params, data, formation_z=None, bin_width=0.5):
    
    z, h, rho_c = params[0], params[1], params[2]
    radial_bins, densities, halo_M_Mean200, halo_R_Mean200 = data[0], data[1], data[2], data[3]
    total_num_halos = halo_M_Mean200.shape[0]
    
    mass_bins = init_mass_bins(halo_M_Mean200, bin_width)
    
    print('')
    if len(mass_bins) != 0:
        # Initialize the results
        num_bins = len(mass_bins) - 1
        results = np.empty((num_bins, 6, args.Nboots))   # [med_mass, Rsp, depth, min_grad, width_dimless, width_physical]
        
        # Bootstrap!
        valid_boots = 0
        reject_times = {ib: 0 for ib in range(num_bins)}
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
            print('Number of halos in bins: ', num_halos_per_bin)
            
            if all(x > 5 for x in num_halos_per_bin):
                
                for i in range(num_bins):
                    if reject_times[i] < 10:
                        # Select halos in the cut and compute density profiles
                        
                        profile_result = fit_density_profile(args, select_radii, select_densities, select_masses, select_r200, 
                                                    10**mass_bins[i], 10**mass_bins[i+1], rho_c)
                        radius, rho, rho_err, fitted_radius, fitted_rho, slope, slope_err, R200_median = profile_result
                        
                        ### If the optimal params are not found! ###
                        if np.all(fitted_rho) == 0:
                            print('This bootstrap is abandoned!')
                            results[:, :, valid_boots] = -1
                            break
                        else:
                            # Fit the slope
                            fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho)) # [dimensionless]
                            
                            # Plot the data and the fit
                            fname = f'form_z{int(np.round(formation_z, 3)*100)}_boots{valid_boots}_bin{i}'
                            if valid_boots == 0:
                                plot_profile(radius, 
                                            rho, rho_err, 
                                            slope, slope_err, 
                                            fitted_radius, fitted_rho, fitted_slope, 
                                            fname, save_dir='result/bootstrap_plots/')
                
                            # Compute the median mass in the mass cut
                            med_mass = compute_median(select_masses, 10**mass_bins[i], 10**mass_bins[i+1])
                            
                            # Compute Rsp
                            physical_fitted_radius = fitted_radius * R200_median
                            Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
                            
                            # depth
                            min_grad = np.min(fitted_slope)
                            min_grad_idx = np.argmin(fitted_slope)
                            print(f'min grad index: {min_grad_idx}')
                            if min_grad_idx < 600:
                                print('This bootstrap is abandoned!')
                                reject_times[i] += 1
                                results[i, :, valid_boots] = -1
                                print('reject times: ', reject_times[i])
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
                    else:
                        print(f'    This mass cut is rejected: {i}')
                        results[i, :, valid_boots] = np.nan
                        continue
            
                ### Only the for loop is complete, update valid_boots
                if np.all(results[:, :, valid_boots] != -1):
                    valid_boots += 1
                    print(f'Nboots updated: {valid_boots}')
                    # reject_times = {i: 0 for i in range(num_bins)}
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
    
    else:
        final_results = {}
        
    return final_results

def plot_profile(radius, 
                 rho, rho_err, 
                 slope, slope_err, 
                 fitted_radius, fitted_rho, fitted_slope, 
                 fname, save_dir=None, save_data=False):
    
    import os
    from matplotlib import pyplot as plt  
    plt.style.use('code/style.mplstyle')
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 15))
    axs[0].scatter(radius, rho, s=1, # color='b', 
                #    label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot$".format(mass_cut[0]+10, mass_cut[1]+10)
                   # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                   )
    axs[0].fill_between(radius, rho-rho_err[:,0], rho+rho_err[:,1], alpha = 0.2, # color = 'b',
                        # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                        )
    axs[0].plot(fitted_radius, fitted_rho, lw=0.5, # color='salmon',
                # label=f'Fit: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                )
    
    # Plot the fitted gradients
    axs[1].scatter(radius, slope, s=1, # color='b',
                  # label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot$".format(mass_cut[0]+10, mass_cut[1]+10)
                   # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                   )
    axs[1].fill_between(radius, slope-slope_err[:,0], slope+slope_err[:,1], alpha = 0.2, # color = 'b',
                        # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                        )
    axs[1].plot(fitted_radius, fitted_slope, lw=0.5, # color='salmon',
                # label=f'Theory: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                )
    
    # General settings
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rho$/$\rho_c$")

    axs[1].set_xscale('log')
    axs[1].set_xlabel(r"r/$R_{200}$")
    axs[1].set_ylabel("Slope")
    axs[1].set_ylim(-6,-0)
    
    # Save the plot
    # start, _ = float_to_str(mass_cut[0], mass_cut[1])'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(save_dir+fname)
    plt.close()
    
    # # Save data
    # if save_data == True:
    #     data = {'radius': radius, 'rho': rho, 'rho_err': rho_err,
    #             'slope': slope, 'slope_err': slope_err, 
    #             'fitted_radius': fitted_radius, 'fitted_rho': fitted_rho, 'fitted_slope': fitted_slope,
    #             'R200_median': R200_median}
        
    #     data_dir = save_dir + f'/data/mass_cut_{start}'
    #     if not os.path.exists(data_dir):
    #        os.makedirs(data_dir)
    #     np.save(os.path.join(data_dir, fname), data)
    
def fit_density_profile(args, radii, densities, masses, r200, bin_start, bin_end, rho_c):
    # Select halos in the cut and compute density profiles
    raw_profiles = stacked_density_profile(radii, densities, masses, r200, 
                                           bin_start, bin_end, rho_c)
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
        print(f'    The profile radius start idx: {start_idx}')
    
    # Compute the slope
    radius, rho, rho_err = filter_profile(radius, rho, rho_err) # Remove zero densities in the centre
    slope = num_deriv(np.log(radius), np.log(rho))                                      # [dimensionless]
    slope_err = num_deriv_err(radius, rho, rho_err)                                     # [dimensionless]
    
    # Fit the density profiles
    fit_profiles = fit_profile_parametric(radius, rho, np.mean(rho_err, axis=1), 1)
    fitted_radius, fitted_rho = fit_profiles[0], fit_profiles[1]
    
    del fit_profiles
    
    return radius, rho, rho_err, fitted_radius, fitted_rho, slope, slope_err, R200_median