""""This script plots depth (top panel) and width (bottom panel) as a function of redshift for a few 
    characteristic mass cuts.
    Hydro simulation only."""

import os
import numpy as np
from func import *
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

print('')
print(f'>>> Plot depth and width vs z <<<')
print('')

root_dir = 'result/bootstrap_stats_phys/'

# ============================================================================================
# Loadd data 
# ============================================================================================

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_list = os.listdir(TNG300_dir)

TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]

TNG300_z, TNG300_bins, TNG300_data = load_data(TNG300_dir, TNG300_snaps)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_list = os.listdir(MTNG_dir)

MTNG_snaps = [264, 237, 214, 179, 151, 129]

MTNG_z, MTNG_bins, MTNG_data = load_data(MTNG_dir, MTNG_snaps)

# Concatenate all data
z = np.concatenate((TNG300_z, MTNG_z))
bins = np.concatenate((TNG300_bins, MTNG_bins))
bins = np.log10(bins)
data = np.concatenate((TNG300_data, MTNG_data), axis=0)

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(2, 1, figsize = (4, 6), dpi=500, sharex=True)
    
uniq_bins = np.unique(bins)
num_bins = len(uniq_bins)
min_bin, max_bin = np.min(bins), np.max(bins)

from matplotlib.colors import LinearSegmentedColormap
red_list = [# '#EE9D9F', '#DE6A69', 
            '#C84747', '#982B2D','#6A0624', '#3D011A'] 
blue_list = ['#89CAEA','#4596CD','#0B75B3','#015696','#012A61','#053061', ]# light to dark
cmap = LinearSegmentedColormap.from_list('my_cmap', blue_list)
# cmap = plt.get_cmap('plasma', num_bins)

bound = np.logspace(min_bin, max_bin+0.1, num_bins) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label('Mass [$M_\\odot$]')
cb._set_scale('log')

# ============================================================================================
# Plot
# ============================================================================================

# Iterate over mass bins
for b in uniq_bins:
    # print('mass cut: ', b)
    idx = np.where(bins == b)
    
    current_z = z[idx]

    current_data = data[idx]
    
    # Plot SORTED data
    sorted_idx = np.argsort(current_z)
    
    X = current_z[sorted_idx]
    
    Y =  current_data[:, 1, 1][sorted_idx]
    Ymin = current_data[:, 1, 1][sorted_idx]-current_data[:, 1, 0][sorted_idx]
    Ymax = current_data[:, 1, 2][sorted_idx]-current_data[:, 1, 1][sorted_idx]
    
    Z = current_data[:, 2, 1][sorted_idx]
    Zmin = current_data[:, 2, 1][sorted_idx]-current_data[:, 2, 0][sorted_idx]
    Zmax = current_data[:, 2, 2][sorted_idx]-current_data[:, 2, 1][sorted_idx]
    
    # Plot depth
    for ix, x in enumerate(X):
        axs[0].plot(x, Y[ix], color=cmap(norm(b)), lw=1, alpha=0.5) 
        axs[0].errorbar(x, Y[ix], yerr=[[Ymin[ix]],[Ymax[ix]]], fmt='.', color=cmap(norm((10**b)/(1+x))))
        # Plot width
        axs[1].plot(x, Z[ix], color=cmap(norm(b)), lw=1, alpha=0.5) 
        axs[1].errorbar(x, Z[ix], yerr=[[Zmin[ix]],[Zmax[ix]]], fmt='.', color=cmap(norm((10**b)/(1+x))))
    
# Final edit
axs[0].set_ylabel("Depth")

axs[1].set_xlabel('z')
axs[1].set_ylabel("Width")

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig5.png')
plt.close()