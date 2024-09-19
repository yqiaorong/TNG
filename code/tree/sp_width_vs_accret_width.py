import numpy as np
import os
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')

DM = '_DM'
snap_list = [8, 13, 25, 40, 67, 99]
acc_width_type = 'percentile' # [std / percentile]


# Load accretion rate width
root_dir = 'result/accretion_rate_plot/'
acc_data = np.load(f'{root_dir}/sim_205_1250{DM}/TNG300{DM}_accret_stats.npy', allow_pickle=True).item()
redshifts = acc_data['redshifts']
# print('redshift: ', np.round(redshifts, 3))
mass_cuts = acc_data['mass_cuts'] # [1.  1.5 2.  2.5 3.  3.5 4. ]
# print('mass cuts: ', mass_cuts)
median_array = acc_data['accret_med']
if acc_width_type == 'std':
    width_array = acc_data['accret_std']
elif acc_width_type == 'percentile':
    lowp_array   = acc_data['accret_low']
    highp_array  = acc_data['accret_high']
    width_array = highp_array - lowp_array
# print(width_array.shape) # (redshifts, mass_cuts) = (6, 6)


# Load splashback features width
Nboots = 1024
feat_idx = 2 # width
root_dir = 'result/bootstrap_stats/'
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250{DM}/Nboots_{Nboots}/'
TNG300_list = os.listdir(TNG300_dir)

# TNG300_z = []
TNG300_all_data = []
for snap in snap_list: # from low z to high z (present)
        
    # Load data
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    z = data['z']
    # TNG300_z.append(z)
    data = data['final_results']
    TNG300_all_data.append(data)
# print('Sanity check of redshifts: ', np.round(TNG300_z, 3))



# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)
if DM == '':
    cmap_name = 'autumn'
elif DM == '_DM':
    cmap_name = 'winter'
cmap = plt.get_cmap(cmap_name, len(TNG300_list))



for cut_idx in range(len(mass_cuts)-1): # from low cut to high cut
    mask = median_array[:, cut_idx] != 0
    plot_x = width_array[mask, cut_idx]
    all_data = [d for d, m in zip(TNG300_all_data, mask) if m==True]

    # Sort all data
    plot_y, plot_y_min, plot_y_max = [], [], []
    
    for snap_idx in range(len(all_data)): 
        if all_data[snap_idx].shape[0] > cut_idx:
            plot_y.append(all_data[snap_idx][cut_idx, feat_idx, 1])
            plot_y_min.append(all_data[snap_idx][cut_idx, feat_idx, 0])
            plot_y_max.append(all_data[snap_idx][cut_idx, feat_idx, 2])
        else:
            pass

    axs.scatter(plot_x[np.argsort(plot_x)], np.array(plot_y)[np.argsort(plot_x)], color=cmap(cut_idx / len(mass_cuts)),
        # label=r'$10^{%.1f}$'%mass_cuts[cut_idx]+'~'
        # +r'$10^{%.1f}$ '%mass_cuts[cut_idx+1]+'$M_\\odot$/h', 
        s=10)
    axs.plot(plot_x[np.argsort(plot_x)], np.array(plot_y)[np.argsort(plot_x)], color=cmap(cut_idx / len(mass_cuts)),
            label=r'$10^{%.1f}$'%mass_cuts[cut_idx]+'~'
            +r'$10^{%.1f}$ '%mass_cuts[cut_idx+1]+'$M_\\odot$/h', 
            linestyle='--')
#     # axs.fill_between(width_array[:len(plot_y), cut_idx], plot_y_min, plot_y_max, 
#     #                 color=cmap(cut_idx / len(mass_cuts)), alpha=0.2)

# Final edit
axs.set_xlabel(f'Accretion rate width ({acc_width_type})')
axs.set_ylabel(f"Splashback feature (width)")
axs.legend()

# Save the plot
save_dir = f'result/width/TNG300{DM}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'TNG300{DM}_sp_width_accret_width_{acc_width_type}'))
plt.close()