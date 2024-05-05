import argparse
import numpy as np
import os
from func import *
from matplotlib import pyplot as plt

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--snap', default=99, type=int)
parser.add_argument('--bin_start', default=0, type=float)
parser.add_argument('--bin_end', default=5.5, type=float)
args = parser.parse_args()

print('')
print(f'>>> Median density profile <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Compute DM softening length
DM_soft = DMsoften(args.snap)
print(f'The DM softening length: {np.round(DM_soft, 3)} ckpc/h')

# Create mass bins
bin_width = 0.5
num_bins = int((args.bin_end-args.bin_start)/bin_width)
mass_bins = np.arange(args.bin_start, args.bin_end+bin_width, bin_width) # mass_bin = x where x: 10^x of 10^10 Msun

# Load directory
load_dir = f'result/DM_halo_density_profiles/snap_{args.snap}'

# Load density profile data
file_list = os.listdir(load_dir)
file_list = [os.path.join(load_dir, fname) for fname in file_list] 

# Set up the final plot
fig, axs = plt.subplots(2, 1, figsize=(10, 15))

# Iteration over mass bins
for i in range(num_bins):
    # Compute median density profiles
    radius, median_rho, rho_err, num_halo = stacked_density_profile(file_list, [mass_bins[i], mass_bins[i+1]], DM_soft)

    if radius.shape != 0:
        # Calculate d log rho / d log r
        slopes = gradient(radius, median_rho)
        slopes_errs = gradient(radius, rho_err)

        # Plot the density profile
        axs[0].errorbar(radius, median_rho, 
                        #yerr = rho_err.T, 
                        fmt='.', 
                    label=f'mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')
        axs[1].errorbar(radius[2:-2], slopes, 
                        #yerr=slopes_errs.T, 
                        fmt='.',
                    label=f'mass bin 10^{mass_bins[i]+10} ~ 10^{mass_bins[i+1]+10} Msun: {num_halo} halos')

axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel("Density [10^10 solar mass / (ckpc/h)^3]")
axs[0].legend()
axs[0].set_title(f'Stacked density profiles')

axs[1].set_xscale('log')
axs[1].set_xlabel("R [r/R_Mean200]")
axs[1].set_ylabel("Slope")
axs[1].legend()
axs[1].set_title('Finding splashback radius')

plt.tight_layout()

# Save directory
save_dir = f'result/Stacked_density_profiles/snap_{args.snap}'
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