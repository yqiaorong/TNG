from matplotlib import pyplot as plt 
plt.style.use('code/plot/style.mplstyle')
import numpy as np
import os
import h5py 
import illustris_python as il

boxsize = 205
res = 1250
data_path = '/n/holylfs05/LABS/hernquist_lab/IllustrisTNG/Runs/'
basePath = data_path + 'L%dn%dTNG/output'%(boxsize,res)



# Save directory
data_path = f'result/bootstrap'
save_dir = f'{data_path}/sim_{boxsize}_{res}'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)



snaps = [33, 40, 50, 67, 78, 99]
mass_cut = [10**12, 10**12.5, 10**13, 10**13.5]

# feature index [Rsp = 0, depth = 1, width = 2]
feats_idx = 2
feats = ['Rsp', 'depth', 'width']

# Plot 1 Rsp vs mass
fig, axs = plt.subplots(1, 1, 
                        # figsize=(10, 8), 
                        dpi=500)
# axs.set_title(f'Splashback radius vs Mass')

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
    data = np.load(data_path+f'/snap_{snap}_Rsp_stats.npy', allow_pickle = True).item()
    data = data ['final_results']
    
    # axs.errorbar(mass_cut, data[1, :, 0], 
    #              yerr = [np.abs(data[0, :, 0]-data[1, :, 0]), 
    #                      np.abs(data[2, :, 0]-data[1, :, 0])],
    #              color=colours[i], label=f'z = {np.round(z[i], 3)}')
    axs.plot(mass_cut, data[:, feats_idx, 1], label=f'z = {np.round(z[i], 3)}')
    axs.fill_between(mass_cut, data[:, feats_idx, 2], data[:, feats_idx, 0], alpha=0.2)
    
axs.set_xlabel('Mass [$M_\\odot$/h]')
axs.set_ylabel(r"$R_{sp}$ [kpc]")
axs.set_xscale('log')
axs.set_yscale('log')
axs.legend()
# plt.tight_layout() # incompatible with pltstyle
plt.savefig(os.path.join(save_dir, f'dm_{feats[feats_idx]}_vs_mass_TNG300'))
plt.close()