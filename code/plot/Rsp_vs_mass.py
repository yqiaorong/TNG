import h5py
import os
import numpy as np
from matplotlib import pyplot as plt
from func import bootstrap

root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250


# Set up the final plot
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
axs.set_title(f'Splashback radius vs Mass')



# Mass bins 
mass_bins = [10**12, 10**12.5, 10**13, 10**13.5]

# Read data
data_path = f'data/{root_dir}/sim_{boxsize}_{res}'
with h5py.File(os.path.join(data_path,'DMhalo_profiles.hdf5'), 'r') as data_f:
    
    snaps = data_f.keys()
    
    # Set up the colour range
    cmap = plt.cm.get_cmap('hsv')
    colours = [cmap(i / len(snaps)) for i in range(len(snaps))]
    
    # Iterate over snapshots
    for s, c in zip(snaps, colours):
        # Import data
        z = np.round(data_f[f'{s}/z'], 3)
        slopes_fit = data_f[f'{s}/slopes_fit']
        radius_fit = data_f[f'{s}/radius_fit']
 
        # Bootstrap the gradients
        num_stacks = slopes_fit.shape[0]
        boots_grads = np.array([bootstrap(slopes_fit[stack_idx], np.min, 100) 
                         for stack_idx in range(num_stacks)])
        # Select 16, 50, 84 percentiles of gradients
        final_grads = np.array([np.percentile(boots_grads[stack_idx], [16, 50, 84]) 
                         for stack_idx in range(num_stacks)])
        # Find Rsp with errors [kpc]
        differ = np.array([np.array([np.abs(slopes_fit[stack_idx]-final_grads[stack_idx, i]) 
                            for stack_idx in range(num_stacks)]) for i in range(3)])
        Rsp_idx = np.array([np.array([np.argmin(differ[i, stack_idx]) 
                            for stack_idx in range(num_stacks)]) for i in range(3)])
        Rsp_with_errs = np.array([np.array([radius_fit[stack_idx, Rsp_idx[i, stack_idx]] 
                                 for stack_idx in range(num_stacks)]) for i in range(3)]).T
        
        # Plot
        axs.errorbar(mass_bins[:num_stacks], Rsp_with_errs[:,1], 
                     yerr = [abs(Rsp_with_errs[:,1]-Rsp_with_errs[:,0]), 
                             abs(Rsp_with_errs[:,2]-Rsp_with_errs[:,1])], 
                     color=c, label=f'z = {z}')
        # axs.plot(mass_bins[:num_stacks], Rsp_with_errs[:,1], color=c, label=f'z = {z}')
        # axs.fill_between(mass_bins[:num_stacks], Rsp_with_errs[:,0], Rsp_with_errs[:,2],
        #                  alpha= 0.5, color=c, label=f'z = {z}')
        
axs.set_xlabel('Mass [MSun/h]')
axs.set_ylabel(r"$R_{sp}$[kpc]")
axs.set_xscale('log')
axs.legend()
plt.tight_layout()
    
# Save directory
save_dir = f'result/Rsp_{root_dir}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'Rsp_vs_mass'))