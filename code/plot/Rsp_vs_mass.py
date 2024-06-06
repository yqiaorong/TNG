import h5py
import os
import numpy as np
from matplotlib import pyplot as plt

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
        # Redshfit
        z = np.round(data_f[f'{s}/z'], 3)
        # Index of Rsp
        indices = np.argmin(data_f[f'{s}/slopes_fit'], axis=1)
        # Rsp
        Rsp = [data_f[f'{s}/radius_fit'][i, idx] for i, idx in enumerate(indices)] # [kpc]
        
        # Plot
        axs.plot(mass_bins[:len(Rsp)], Rsp, color=c, label=f'z = {z}')
        
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