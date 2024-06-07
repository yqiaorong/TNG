import h5py
import os
import numpy as np
from matplotlib import pyplot as plt

root_dir = 'DMhalo_density_profiles_old'
boxsize = 205
res = 1250



# Set up the final plot
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
axs.set_title(f'Splashback radius vs Redshift')



# Mass bins 
mass_bins = [12, 12.5, 13, 13.5]

# Set up the colour range
cmap = plt.cm.get_cmap('hsv')
colours = [cmap(i / len(mass_bins)) for i in range(len(mass_bins))]

# Read data
data_path = f'data/{root_dir}/sim_{boxsize}_{res}'
with h5py.File(os.path.join(data_path,'DMhalo_profiles.hdf5'), 'r') as data_f:
    
    snaps = data_f.keys()
    
    # # Set up the colour range
    # cmap = plt.cm.get_cmap('hsv')
    # colours = [cmap(i / len(snaps)) for i in range(len(snaps))]
    
    # Iterate over snapshots
    Rsp_all, z_all = [], []
    for s in snaps:
        # Redshfit
        z = np.round(data_f[f'{s}/z'], 3)
        z_all.append(z)
        
        # Index of Rsp
        indices = np.argmin(data_f[f'{s}/slopes_fit'], axis=1)
        # Rsp
        Rsp = [data_f[f'{s}/radius_fit'][i, idx] for i, idx in enumerate(indices)] # [kpc]
        Rsp_all.append(Rsp)

for i in range(len(mass_bins)):
    Rsp_bin = []
    for sublist in Rsp_all:
        if len(sublist) >= i+1:
            Rsp_bin.append(sublist[i])
    
    # Plot
    axs.plot(z_all[-len(Rsp_bin):], Rsp_bin, color=colours[i], label=f'mass = 10^{mass_bins[i]} MSun/h')
        
axs.set_xlabel('z')
axs.set_ylabel(r"$R_{sp}$[kpc]")
axs.legend()
plt.tight_layout()
    
# Save directory
save_dir = f'result/Rsp_{root_dir}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'Rsp_vs_redshift'))