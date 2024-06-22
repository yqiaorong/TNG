import argparse
import numpy as np
import os
from func import *
from matplotlib import pyplot as plt  

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',   default=205,  type=int)
parser.add_argument('--res',       default=1250, type=int)
parser.add_argument('--snapnum',   default=99,   type=int)
parser.add_argument('--bin_start', default=3.5,  type=float)
parser.add_argument('--bin_end',   default=4,    type=float)

parser.add_argument('--root_dir',  default=None, type=str)
args = parser.parse_args()

print('')
print(f'>>> Median density profile <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')


# Create mass bins
bin_width = 0.5
num_bins = int((args.bin_end-args.bin_start)/bin_width)
mass_bins = np.arange(args.bin_start, args.bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun/h

# Load directory
load_dir = 'result/'+args.root_dir+f'/sim_{args.boxsize}_{args.res}/snap_{args.snapnum}/densities'

# Load density profile data
file_list = os.listdir(load_dir)
file_list = [os.path.join(load_dir, fname) for fname in file_list] 

# Set up the final plot
fig, axs = plt.subplots(2, 1, figsize=(10, 15))
        
# Iteration over mass bins
for i in range(num_bins):
    
    # Compute median density profiles
    raw_profiles = stacked_density_profile(file_list, [mass_bins[i], mass_bins[i+1]])
    radius, rho, rho_err = raw_profiles[0][1:], raw_profiles[1][1:], raw_profiles[2][1:] # [dimensionless]
    num_halo, R200_median = raw_profiles[3], raw_profiles[4]                    # [ckpc/h]
    del raw_profiles
    
    # Cross check if rho is zero
    mask = rho != 0
    radius, rho, rho_err = radius[mask], rho[mask], rho_err[mask]
    
    # Calculate d log rho / d log r
    # grad_result = gradient(radius, median_rho, rho_err[:,0])
    # slopes_radius, slopes, slopes_errs = grad_result[0], grad_result[1], grad_result[2] 
    slopes = num_deriv(np.log(radius), np.log(rho))              # [dimensionless]
    slopes_err = num_deriv_err(radius, rho, rho_err)


    
    # Fit the density profiles
    fit_profiles = fit_profile_parametric(radius, rho, np.mean(rho_err, axis=1), 1)
    new_radius, new_rho = fit_profiles[0], fit_profiles[1]
    del fit_profiles
    
    # Fit the slope 
    # M1
    # fit_results = fit_gradient_parametric(radius[20:], rho[20:], np.mean(rho_err, axis=1)[20:],
    #                                       slopes[20:], np.mean(slopes_err, axis=1)[20:], 1)
    # slopes_fit_r, slopes_fit = fit_results[0], fit_results[1]
    # del fit_results
    # M2
    slopes_fit = num_deriv(np.log(new_radius), np.log(new_rho))      # [dimensionless]
    slopes_fit_r = new_radius                                        # [dimensionless]
    
    

    # Plot the density profile
    axs[0].scatter(radius, rho, s=1, color='b',
                   label=f'Data: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    axs[0].fill_between(radius, rho-rho_err[:,0], rho+rho_err[:,1], alpha = 0.2, color = 'b',
                        label=f'Errorbar: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    axs[0].plot(new_radius, new_rho, lw=0.5, color='salmon',
                label=f'Fit: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    
    # Plot the fitted gradients
    axs[1].scatter(radius, slopes, s=1, color='b',
                   label=f'Data: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    axs[1].fill_between(radius, slopes-slopes_err[:,0], slopes+slopes_err[:,1], alpha = 0.2, color = 'b',
                        label=f'Errorbar: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    axs[1].plot(slopes_fit_r, slopes_fit, lw=0.5, color='salmon',
                label=f'Theory: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun/h: {num_halo} halos')
    
    

    # General settings
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rho$/$\rho_c$")
    axs[0].legend()
    axs[0].set_title(f'Stacked density profiles')

    axs[1].set_xscale('log')
    axs[1].set_xlabel("r/R200")
    axs[1].set_ylabel("Slope")
    axs[1].set_ylim(-6,-0)
    axs[1].legend()
    axs[1].set_title('Finding splashback radius')

plt.tight_layout()



# Save directory
save_dir = f'result/Stacked_{args.root_dir}/sim_{args.boxsize}_{args.res}/snap_{args.snapnum}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Plot name and save
# start, end = float_to_str(args.bin_start, args.bin_end)
start, end = float_to_int(args.bin_start, args.bin_end)
plt.savefig(os.path.join(save_dir, f'Bins_{start}_to_{end}'))

# Save data
save_data = {'profile_radius': radius,           'profile_densities': rho,           # [dimensionless]
             'profile_radius_fit': new_radius,   'profile_densities_fit': new_rho,   # [dimensionless]
             'slopes_radius': radius,            'slopes': slopes,                   # [dimensionless]
             'slopes_radius_fit': slopes_fit_r,  'slopes_fit': slopes_fit,           # [dimensionless]
             'R200_median': R200_median, # [ckpc/h]
             'density_err': rho_err,             'slope_err': slopes_err             # [dimensionless]
             }                                             
np.save(os.path.join(save_dir, f'Bins_{start}_to_{end}'), save_data)