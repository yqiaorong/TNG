# If func has been imported, ignore the import
import numpy as np
import math
import os
from func import *
from func2 import *

def init_bins(args, array_1d, bin_type, bins=None, bin_width=None, bin_start=None, bin_end=None):
    
    # Set up the bins
    if bins is None:
        if bin_type == 'mass':
            # Logspace
            bin_start, bin_end, bin_width = 3, 6, 0.5
            bins = np.logspace(bin_start, bin_end, num=int((bin_end - bin_start) / bin_width) + 1)
        else:
            if bin_start is None or bin_end is None:
                min_bin, max_bin = np.min(array_1d), np.max(array_1d)
                bin_start, bin_end = math.floor(min_bin*2)/2, math.ceil(max_bin*2)/2
            if bin_width is None:
                # Ask to enter bin_width
                bin_width = float(input(f'Enter bin width for {args.bin_type}: '))
            bins = np.arange(bin_start, bin_end+bin_width, bin_width)
            print('old bins: ', bins)

    # Check the number of halos in the largest bin. if it's less than 2, reject this bin
    new_bins = []
    for ib, b in enumerate(bins[:-1]):
        num_halos_in_bin = np.where((array_1d >= b) & (array_1d < bins[ib+1]))[0].shape[0]
        if num_halos_in_bin > 50:
            new_bins.append(b)
            
    if len(new_bins) != 0:
        if args.bin_type == 'mass':
            new_bins.append(10**(np.log10(new_bins[-1])+bin_width))
        else:
            new_bins.append(new_bins[-1]+bin_width)
    else:
        print('No valid bins found!')
        exit()
    del bins
    
    print('new bins: ', new_bins)
    return new_bins
    
    
# Bootstrap setup
def bootstrap_all_features(args, params, data, x_type, 
                           plot_x_type_name=None, 
                           x_bins=None, x_width=None, x_start=None, x_end=None):
    
    # Plot the data and the fit
    if plot_x_type_name is None:
        plot_x_type_name = x_type
    
    # Initialise
    z, h = params[0], params[1]
    radial_bins, densities, x_vals, halo_R_Mean200 = data[0], data[1], data[2], data[3]
    total_num_halos = x_vals.shape[0]

    x_bins = init_bins(args, x_vals, bin_type=x_type, bins=x_bins, bin_width=x_width, bin_start=x_start, bin_end=x_end)
    
    if len(x_bins) != 0:
        # Initialize the results
        num_bins = len(x_bins) - 1
        results = np.empty((num_bins, 7, args.Nboots))   # [med_bin_type, Rsp, depth, min_grad, width_dimless, width_physical, DWratio]
        if args.profile == 'DK14':
            params_results = np.empty((num_bins, 10, args.Nboots)) 
        elif args.profile == 'functional':
            params_results = np.empty((num_bins, 5, args.Nboots)) 
        # Bootstrap!
        valid_boots = 0
        reject_times = {ib: 0 for ib in range(num_bins)}
        while valid_boots < args.Nboots: 
            
            # Random selection of halos
            indices = np.random.randint(0, total_num_halos, args.Nsample)
            # Select densities and masses
            select_radii     = radial_bins[indices]    # [kpc]
            select_densities = densities[indices]      # [Msun / (kpc)^3]
            select_vals      = x_vals[indices]       
            select_r200      = halo_R_Mean200[indices] # [kpc]
            del indices

            # Calculating the number of halos in each cut
            num_halos_per_bin = count_halos(select_vals, x_bins[:-1], x_bins[1:])
            print(f'Number of halos in {x_type} bins: ', num_halos_per_bin)

            if all(x > 5 for x in num_halos_per_bin):
                
                for i in range(num_bins):
                    print(f'>>>>>> Processing {x_type} cut {i}: <<<<<<')
                    if reject_times[i] < args.reject_limit:
                        # Select halos in the cut and compute density profiles
                        
                        fit_data = fit_density_profile(args, select_radii, select_densities, select_vals, select_r200, 
                                                             x_bins[i], x_bins[i+1])
                        fit_profile, fit_slope, R200_median, fitted_params = fit_data[0], fit_data[1], fit_data[2], fit_data[3] 
                        radius, rho, rho_err, slope, slope_err = fit_profile[0], fit_profile[1], fit_profile[2], fit_profile[3], fit_profile[4]
                        fitted_radius, fitted_rho, fitted_slope = fit_slope[0], fit_slope[1], fit_slope[2]
                        del fit_data, fit_profile, fit_slope
                        # Fit the slope
                        fitted_slope = num_deriv(np.log(fitted_radius), np.log(fitted_rho)) # [dimensionless]
                        
                        # If the optimal params are not found! 
                        if np.all(fitted_slope) == 0:
                            print('     This bootstrap is abandoned!-----No fitted profile found!')
                            results[:, :, valid_boots] = -1
                            print('reject times: ', reject_times[i])
                            break
                                        
                        # If the left part of the gradients are not increasing, keep this bootstrap
                        # if not np.all(np.diff(fitted_slope[fitted_radius < 0.1]) < 1E-3):
                        #     print('     This bootstrap is abandoned!-----Left part of the gradients are increasing!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     # exit()
                        #     break                          
        
                        # Compute the median mass in the mass cut
                        med = compute_median(select_vals, x_bins[i], x_bins[i+1])
                        
                        # Compute Rsp
                        physical_fitted_radius = fitted_radius * R200_median
                        Rsp_dimless = fitted_radius[np.argmin(fitted_slope)]  # [dimensionless]
                        Rsp = physical_fitted_radius[np.argmin(fitted_slope)] # [kpc]
                        
                        # # If the gradients have more than one minimum, reject this bootstrap
                        # if test_shape(fitted_radius, fitted_slope, Rsp_dimless):
                        #     print('     This bootstrap is abandoned!-----More than one minimum found!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     break
                            
                        # depth
                        min_grad = np.min(fitted_slope)
                        min_grad_idx = np.argmin(fitted_slope)
                        print(f'min grad index: {min_grad_idx}, min grad val: {min_grad}')
                        # if args.bin_type == 'NFWconc' and i == 3 and min_grad < -4:
                        #     print('     This bootstrap is abandoned!-----Inaccurate Rsp found!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     break
                        # elif args.bin_type == 'mass' and min_grad < -5:
                        #     print('     This bootstrap is abandoned!-----Inaccurate Rsp found!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     break
                        
                        # if min_grad_idx < 600:
                        #     print('     This bootstrap is abandoned!-----Inaccurate Rsp found!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     break
                   
                        left_data = fitted_slope[:min_grad_idx]
                        right_data = fitted_slope[min_grad_idx:] 
                        
                        max_grad = np.max(right_data) 
                        if max_grad > 0:
                            max_grad = np.max(slope) 
                        if max_grad > 0:
                            max_grad = 0
                        print('max grad', max_grad)
                        depth = max_grad - min_grad
                        
                        # Width
                        half_grad = min_grad + depth/2
                        
                        left_idx = np.argmin(np.abs(left_data - half_grad))
                        print(f'left width index: {left_idx}')
                        
                        right_idx = min_grad_idx + np.argmin(np.abs(right_data - half_grad))
                        
                        width = physical_fitted_radius[right_idx] - physical_fitted_radius[left_idx]
                        width_dimless = fitted_radius[right_idx] - fitted_radius[left_idx]
                        
                        # Depth vs width
                        DWratio = depth / width_dimless
                        # print(f'DW ratio: {DWratio}')
                        # if DWratio > 8 or DWratio < 0.5:
                        #     print('     This bootstrap is abandoned!-----Inaccurate DWratio found!')
                        #     reject_times[i] += 1
                        #     results[i, :, valid_boots] = -1
                        #     print(f'     Cut {i} reject times: ', reject_times[i])
                        #     break
                        
                        try:
                            plot_dir = f'result/bootstrap_plots_{args.profile}/{args.sim}/with_{plot_x_type_name}/snap_{args.snapnum}/cut_{i}/'
                        except AttributeError:
                            plot_dir = f'result/bootstrap_plots_{args.profile}/with_{plot_x_type_name}/cut_{i}/'
                        if not os.path.exists(plot_dir):
                            os.makedirs(plot_dir)
                        plot_fname = f'boot_{valid_boots}.png'
                        plot_profile(radius=radius, rho=rho, rho_err=rho_err, 
                                     slope=slope, slope_err=slope_err, 
                                     fitted_radius=fitted_radius, fitted_rho=fitted_rho, fitted_slope=fitted_slope,
                                        depth_coords = [(Rsp_dimless, min_grad), (Rsp_dimless, max_grad)],
                                        width_coords = [(fitted_radius[left_idx], half_grad), (fitted_radius[right_idx], half_grad)], 
                                        plot_dir=plot_dir, plot_fname=plot_fname, plot_text=f'cut_{i}_{plot_fname.split('_')[1]}')
            
                        # Append results
                        results[i, :, valid_boots]        = med, Rsp, depth, min_grad, width_dimless, width, DWratio
                        params_results[i, :, valid_boots] = fitted_params
                        
                    else:
                        print(f'>>>>>>>>>> This {x_type} cut {i} is rejected: {reject_times[i]} <<<<<<<<<<')
                        results[i, :, valid_boots]        = np.nan
                        params_results[i, :, valid_boots] = np.nan
            
                ### Only the for loop is complete, update valid_boots
                if np.all(results[:, :, valid_boots] != -1):
                    valid_boots += 1
                    print(f'Nboots updated: {valid_boots}')
                else:
                    print(f'Nboots not updated: {valid_boots} ')
                print('')
                
        # Get the statistical results 
        final_results = {'z': z, 'h': h, f'{x_type}_bins': x_bins[:-1]}
        final_results[f'med_{x_type}'] = np.percentile(results[:, 0, :], [16, 50, 84], axis=1)
        final_results['Rsp']            = np.percentile(results[:, 1, :], [16, 50, 84], axis=1)
        final_results['depth']          = np.percentile(results[:, 2, :], [16, 50, 84], axis=1)
        final_results['abs_depth']      = np.percentile(results[:, 3, :], [16, 50, 84], axis=1)    
        final_results['width_dimless']  = np.percentile(results[:, 4, :], [16, 50, 84], axis=1)
        final_results['width_physical'] = np.percentile(results[:, 5, :], [16, 50, 84], axis=1) 
        final_results['DWratio']        = np.percentile(results[:, 6, :], [16, 50, 84], axis=1)
        final_results['full_results']   = results 
        print(final_results.keys())
        
        # Final params results
        final_params_results = {}
        if args.profile == 'DK14':
            final_params_keys = ['rho_s', 'r_s', 'alpha', 'r_t', 'beta', 'gamma', 'rho_g', 'b_e', 'S_e', 'R200']
        elif args.profile == 'functional':
            final_params_keys = ['rho_s', 'r_s', 'alpha', 'r_t', 'beta']
        for i, key in enumerate(final_params_keys):
            final_params_results[key] = np.percentile(params_results[:, i, :], [16, 50, 84], axis=1)
        final_results['params'] = final_params_results
        print(final_params_results.keys())
    
    else:
        final_results = {}
        final_params_results = {}
        
    return final_results, final_params_results


def test_shape(radius, slopes, Rsp):
    from scipy.ndimage import median_filter
    
    ctr_idx = np.where(radius == Rsp)[0][0]
    middles = slopes[int(ctr_idx-300) : int(ctr_idx+400)]
    smoothed = median_filter(middles, size=11)
    residual = middles - smoothed

    if np.any(np.abs(residual) > 0.05):
        return True
    else:
        return False


def plot_profile(radius, rho, slope, 
                 fitted_radius, fitted_rho, fitted_slope, 
                 plot_dir, plot_fname, plot_text=None,
                 rho_err=None, slope_err=None, 
                 depth_coords=None, width_coords=None,):
    
    import os
    from matplotlib import pyplot as plt  
    plt.style.use('code/style.mplstyle')
    
    fig, axs = plt.subplots(2, 1, figsize=(2, 3))
 
    axs[0].scatter(radius, rho, s=1, # color='b', 
                #    label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot$".format(mass_cut[0]+10, mass_cut[1]+10)
                   # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                   )
    if rho_err is not None:
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
    if slope_err is not None:
        axs[1].fill_between(radius, slope-slope_err[:,0], slope+slope_err[:,1], alpha = 0.2, # color = 'b',
                            # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                            )
    axs[1].plot(fitted_radius, fitted_slope, lw=0.5, # color='salmon',
                # label=f'Theory: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                )
    
    # Plot the depth and width as segments on axs1
    if depth_coords is not None:
        x_depth, y_depth = zip(*depth_coords)
        x_width, y_width = zip(*width_coords)
        axs[1].plot(x_depth, y_depth, color='red', linestyle='-', label='Depth')
        axs[1].plot(x_width, y_width, color='blue', linestyle='-', label='Width')
    
    # Plot text
    if plot_text is not None:
        axs[0].text(0.5, 0.8, plot_text, transform=axs[0].transAxes, fontsize=8, verticalalignment='top')
    
    # General settings
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rho$/$\rho_c$")

    axs[1].set_xscale('log')
    axs[1].set_xlabel(r"r/$R_{200}$")
    axs[1].set_ylabel("Slope")
    axs[1].set_ylim(-6,-0)
    
    # Save the plot
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    plt.savefig(plot_dir+plot_fname, dpi=100)
    plt.close()
    

def fit_density_profile(args, radii, densities, targets, r200, bin_start, bin_end):

    # Select halos in the cut and compute density profiles
    raw_profiles = stacked_density_profile(radii, densities, targets, r200, 
                                           bin_start, bin_end)
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
        print(f'The profile radius start idx: {start_idx}')
        
    # Compute the slope
    radius, rho, rho_err = filter_profile(radius, rho, rho_err) # Remove zero densities in the centre
    slope = num_deriv(np.log(radius), np.log(rho))              # [dimensionless]
    slope_err = num_deriv_err(radius, rho, rho_err)             # [dimensionless]
    
    
    if args.profile == 'DK14':
        fit_profile = fit_profile_parametric_DK14
        fit_gradient = fit_gradient_parametric_DK14
    elif args.profile == 'functional':
        fit_profile = fit_profile_parametric_fk
        fit_gradient = fit_gradient_parametric_fk
    
    # Fit the density profiles
    fit_profiles = fit_profile(radius, rho, 1, # density_errors=np.mean(rho_err, axis=1)
                               )
    fitted_radius, fitted_rho, fitted_params = fit_profiles[0], fit_profiles[1], fit_profiles[2] 
    del fit_profiles
    
    # Fit the slope
    fit_slopes = fit_gradient(radius, slope, 1, init_p0=fitted_params, # gradients_errors=np.mean(slope_err, axis=1)
                              )
    fitted_slope, fitted_params = fit_slopes[1], fit_slopes[2]
    del fit_slopes
    
    return [radius, rho, rho_err, slope, slope_err], [fitted_radius, fitted_rho, fitted_slope], R200_median, fitted_params