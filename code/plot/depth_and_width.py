import h5py
import os
from matplotlib import pyplot as plt
import numpy as np

root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250



# Set up the final plot
fig, axs = plt.subplots(2, 1, figsize=(10, 12))
axs[0].set_title(f'Splashback feature depth')
axs[1].set_title(f'Splashback feature width')



# Mass bins 
mass_bins = [10**11.5, 10**12, 10**12.5, 10**13, 10**13.5]



data_path = f'data/{root_dir}/sim_{boxsize}_{res}'
with h5py.File(os.path.join(data_path,'DMhalo_profiles.hdf5'), 'r') as data_f:
    
    snaps = list(data_f.keys())
    
    # Set up the colour range
    cmap = plt.cm.get_cmap('hsv')
    colours = [cmap(i / len(snaps)) for i in range(len(snaps))]
    
    # Iterate over snapshots
    for s, c in zip(snaps[:6], colours[:6]):
        
        # Import data
        z = np.round(data_f[f'{s}/z'], 3)
        slopes_fit = data_f[f'{s}/slopes_fit']
        radius_fit = data_f[f'{s}/radius_fit']

        depths, widths = [], []
        num_stacks = slopes_fit.shape[0]
        for row in range(num_stacks):
            
            # Depth
            min_grad = np.min(slopes_fit[row])
            min_grad_idx = np.argmin(slopes_fit[row])
            
            left_data = slopes_fit[row, :min_grad_idx]
            right_data = slopes_fit[row, min_grad_idx:]
            
            max_grad = np.max(right_data)
            depth = max_grad - min_grad
            depths.append(depth)
            
            # Width
            half_grad = min_grad + depth/2
            
            left_idx = np.argmin(np.abs(left_data - half_grad))
            right_idx = min_grad_idx + np.argmin(np.abs(right_data - half_grad))
            
            width = radius_fit[row, right_idx] - radius_fit[row, left_idx]
            widths.append(width)

        # Plot 
        axs[0].plot(mass_bins[:num_stacks], depths, color=c, label=f'z = {z}')
        axs[1].plot(mass_bins[:num_stacks], widths, color=c, label=f'z = {z}')
        
axs[1].set_xlabel('Mass [MSun/h]')
axs[0].set_ylabel('Depth')
axs[1].set_ylabel('Width [kpc]')
axs[0].set_xscale('log')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[0].legend()
axs[1].legend()
plt.tight_layout()

# Save directory
save_dir = f'result/Rsp_{root_dir}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'Depth and Width'))
        
        