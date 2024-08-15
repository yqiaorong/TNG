from matplotlib import pyplot as plt
plt.style.use('code/style.mplstyle') 
import numpy as np
import os
import h5py 
import illustris_python as il

boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Save directory
data_path = f'result/bootstrap_stats/sim_{boxsize}_{res}'
save_dir = f'result/bootstrap_plot/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)



snaps = [33, 40, 50, 67, 78, 99]
mass_cut = [10**12, 10**12.5, 10**13, 10**13.5]

# Plot 1 Rsp vs mass
fig, axs = plt.subplots(2, 1, dpi=400,
                        # figsize=(10, 12)
                        )
axs[0].set_title(f'Splashback features (Depth) vs Mass')
axs[1].set_title(f'Splashback features (Width) vs Mass')

# Set up the colour range
# cmap = plt.cm.get_cmap('hsv')
# colours = [cmap(i / len(snaps)) for i in range(len(snaps))]

z = []
for i, snap in enumerate(snaps):
    # Load redshifts
    with h5py.File(il.snapshot.snapPath(basePath, snap), 'r') as f:
        header = dict(f['Header'].attrs.items())
        scale_factor = header['Time']
        z.append(1 / scale_factor - 1)
    # Load data
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    data = data ['final_results']
    
    depth_idx, width_idx = 1, 2
    # Plot 1: depth
    # M1
    axs[0].errorbar(mass_cut, data[:, depth_idx, 1], 
                 yerr = [np.abs(data[:, depth_idx, 0]-data[:, depth_idx, 1]),  # Lower bound
                         np.abs(data[:, depth_idx, 2]-data[:, depth_idx, 1])], # Upper bound 
                 label=f'z = {np.round(z[i], 3)}')
    # M2
    # axs[0].plot(mass_cut, data[:, depth_idx, 1], label=f'z = {np.round(z[i], 3)}')
    # axs[0].fill_between(mass_cut, data[:, depth_idx, 2], data[:, depth_idx, 0], alpha=0.2,  
    #                  # label=f'z = {np.round(z[i], 3)}'
    #                  )

    # Plot 2: width
    # M1
    axs[1].errorbar(mass_cut, data[:, width_idx, 1], 
                 yerr = [np.abs(data[:, width_idx, 0]-data[:, width_idx, 1]),  # Lower bound
                         np.abs(data[:, width_idx, 0]-data[:, width_idx, 2])], # Upper bound
                 label=f'z = {np.round(z[i], 3)}') 
    # M2 
    # axs[1].plot(mass_cut, data[:, width_idx, 1], label=f'z = {np.round(z[i], 3)}')
    # axs[1].fill_between(mass_cut, data[:, width_idx, 0], data[:, width_idx, 2], alpha=0.2,  
    #                  # label=f'z = {np.round(z[i], 3)}'
    #                  )
    
axs[1].set_xlabel('Mass [$M_\\odot$/h]')
axs[0].set_ylabel('Depth')
axs[1].set_ylabel('Width')

axs[0].set_xscale('log')
axs[1].set_xscale('log')

axs[0].legend(loc='best')
axs[1].legend(loc='best')
plt.savefig(os.path.join(save_dir, f'Rsp_features'))
plt.close()