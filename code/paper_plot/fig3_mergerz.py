import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, get_bins, plot_feature

snap = 264
bin_type = 'mergerz'
xlabel = r'$z_{\rm merger}$'
    
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

simus = 'Hydro'
features = ['depth', 'width_dimless', 'DWratio']                                                                                                                     
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$", r"$\mathcal{D}/\mathcal{W}$"]

# ============================================================================================
# Set up the plot
# ============================================================================================

import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(9, 5), dpi=500)
gs = gridspec.GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 0.05], hspace=0.05)  # 2 rows for plots, 1 for colorbar

axs = np.empty((2, 3), dtype=object)  # 2 rows × 3 columns subplot array

# Create subplots for top 2 rows
for col in range(3):
    # axs[row, col] = fig.add_subplot(gs[row, col])
    axs[0, col] = fig.add_subplot(gs[0, col])
    axs[1, col] = fig.add_subplot(gs[1, col], sharex=axs[0, col]) 
    axs[0, col].tick_params(labelbottom=False)

# ============================================================================================
# Set up the colorbar
# ============================================================================================

# Use mass cuts
min_bin, max_bin, bin_width = 13, 15, 0.25
num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

cmap = plt.get_cmap('plasma', num_bins)
bound = np.logspace(min_bin, max_bin+0.01, num_bins+1) 
norm = BoundaryNorm(bound, cmap.N)

# Create long horizontal colorbar in the third row
cbar_ax = fig.add_subplot(gs[2, :])  # full-width colorbar
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax,
                  orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label(r'$M_{200m}/M_{\odot}$')
cb._set_scale('log')

# ============================================================================================
# Plot
# ============================================================================================

for ifeat, feature in enumerate(features):
    # # Plot colour bar
    # cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs[1,ifeat], 
    #                   orientation='horizontal', spacing='proportional', ticks=bound)
    # cb.set_label(r'$M_{200m}/M_{\odot}$')
    # cb._set_scale('log')
   
    # Row 0
    MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(f'result/bootstrap_stats_DK14/with_{bin_type}/MTNG/{simus}-Arepo/MTNG-L500-4320-A/Nboots_1024/',
                                                      f'snap_{snap}_Rsp_stats.npy',
                                                    f'med_{bin_type}', feature)            
    plot_feature(simus, MTNG_z, MTNG_bin_data, MTNG_feat, [axs[0,ifeat], cmap, norm], label=f'z={MTNG_z:.1f}')
    
    # Row 1
    data_dir = f'result/bootstrap_stats_DK14/with_{bin_type}_perMassCut/MTNG/{simus}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
    for fname in os.listdir(data_dir):
        print(fname)
        _, MTNG_bin_data, MTNG_feat = load_stats(data_dir, fname, f'med_{bin_type}', feature)    
        bin_val = 10**(10+float(fname.split('_')[3])/10)
        plot_feature(simus, bin_val, MTNG_bin_data, MTNG_feat, [axs[1,ifeat], cmap, norm])
    
    axs[0,ifeat].set_ylabel(Ylabels[ifeat])
    axs[1,ifeat].set_ylabel(Ylabels[ifeat])
    axs[1,ifeat].set_xlabel(xlabel)

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/fig3/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig3_{bin_type}_{simus}')
plt.close()