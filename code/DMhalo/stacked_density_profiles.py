import argparse
import numpy as np
import os
from func import *
from matplotlib import pyplot as plt

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--boxsize',default=205,type=int)
parser.add_argument('--res',default=1250,type=int)
parser.add_argument('--snapnum', default=99, type=int)
parser.add_argument('--bin_start', default=3.5, type=float)
parser.add_argument('--bin_end', default=4, type=float)
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
mass_bins = np.arange(args.bin_start, args.bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun

# Load directory
load_dir = f'result/DMhalo_density_profiles/sim_{args.boxsize}_{args.res}/snap_{args.snapnum}/densities'

# Load density profile data
file_list = os.listdir(load_dir)
file_list = [os.path.join(load_dir, fname) for fname in file_list] 

# Set up the final plot
fig, axs = plt.subplots(2, 1, figsize=(10, 15))
        
# Iteration over mass bins
for i in range(num_bins):
    # Compute median density profiles
    radius, median_rho, rho_err, num_halo = stacked_density_profile(file_list, [mass_bins[i], mass_bins[i+1]])
    # Compute fitted median density profiles
    result = fit_profile_parametric(radius, median_rho, rho_err[:,0], 1)
    
    # Calculate d log rho / d log r
    slope_radius, slopes, slopes_errs = gradient(radius, median_rho, rho_err[:,0])
    
    # Fit the slope
    # M1:
    slopes_fits_r, slopes_fits, _ = gradient(result[0], result[1])
    # M2:
    # slope_rhos = median_rho[2:-2]
    # slope_rhos_errs = rho_err[2:-2, 0]
    # slope_result = fit_gradient_parametric(slope_radius, slope_rhos, slope_rhos_errs, slopes, slopes_errs, 1)
    # slopes_fits_r, slopes_fits = slope_result[0], slope_result[1]
    
    # Plot the density profile
    axs[0].errorbar(radius, median_rho, 
                    #yerr = rho_err.T, 
                    fmt='.', 
                label=f'Data: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')
    axs[0].errorbar(result[0], result[1], 
                    #yerr = rho_err.T, 
                    fmt='.', 
                label=f'Fit: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')
    
    # Plot the fitted gradients
    axs[1].errorbar(slope_radius, slopes, 
                        #yerr=slopes_errs.T, 
                        fmt='.',
                    label=f'Data: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')
    axs[1].errorbar(slopes_fits_r, slopes_fits, 
                        #yerr=slopes_errs.T, 
                        fmt='.',
                    label=f'Theory: mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')
      
    # General settings
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylabel("density [M$_{\odot}$/kpc$^3$]")
    axs[0].legend()
    axs[0].set_title(f'Stacked density profiles')

    axs[1].set_xscale('log')
    axs[1].set_xlabel("r/R200")
    axs[1].set_ylabel("Slope")
    axs[1].legend()
    axs[1].set_title('Finding splashback radius')

plt.tight_layout()

# Save directory
save_dir = f'result/Stacked_density_profiles/sim_{args.boxsize}_{args.res}/snap_{args.snapnum}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Plot name and save
if str(args.bin_start).endswith('0'):
    start = int(args.bin_start)
    enda, endb = str(args.bin_end).split('.')
    end = f'{enda}-{endb}'
else:
    starta, startb = str(args.bin_start).split('.')
    start = f'{starta}-{startb}'
    end = int(args.bin_end)
plt.savefig(os.path.join(save_dir, f'Bins_{start}_to_{end}'))