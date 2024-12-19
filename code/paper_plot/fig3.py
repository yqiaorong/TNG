""""This script plots depth (top panel) and width (bottom panel) as a function of mass for a few redshifts.
    Both Hydro and DM-only simulations."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot Depth and width vs mass for a few redshifts <<<')
print('')

root_dir = 'result/bootstrap_stats_phys/'

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(2, 1, figsize=(4, 6), dpi=500, sharex=True)

# ============================================================================================
# Set up the colorbar
# ============================================================================================

TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8
                ]
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

MTNG_snaps = [264, 237, 214, 179, 151, 129]
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)

# Set up the colorbar
# -----------------------------------------------------------------------------------------
from matplotlib.colors import LinearSegmentedColormap
red_list = [# '#EE9D9F', 
            '#DE6A69', '#C84747', '#982B2D','#6A0624', '#3D011A'] # light to dark
blue_list = ['#89CAEA','#4596CD', '#0B75B3', '#015696', '#012A61']
cmap = LinearSegmentedColormap.from_list('my_cmap', red_list)

# cmap = plt.get_cmap('autumn', len(TNG300_snaps))
bound = np.linspace(z_f, z_i, len(TNG300_snaps))
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# -----------------------------------------------------------------------------------------

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.set_label('z')

# ============================================================================================
# Plot
# ============================================================================================

def load_stats(dir, snap):
        data = np.load(f'{dir}/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()

        z = np.round(data['z'], 3)
        mass_bins     = [10**(10+mass) for mass in data['mass_bins']]
        final_results = data['final_results']
        
        return z, mass_bins, final_results

def plot_depth(z, mass_bins, final_results, simu):
        if simu == 'DM':
                ls = '--'
        else:
                ls = '-'
                       
        axs[0].plot(mass_bins, final_results[:, 1, 1], color=cmap(norm(z)), lw=1, alpha=0.5, linestyle=ls) 
        axs[0].errorbar(mass_bins, final_results[:, 1, 1],
                        yerr=[final_results[:, 1, 1]-final_results[:, 1, 0], 
                              final_results[:, 1, 2]-final_results[:, 1, 1]],
                        color=cmap(norm(z)), fmt='.')
        
def plot_width(z, mass_bins, final_results, simu):
        if simu == 'DM':
                ls = '--'
        else:
                ls = '-'
                       
        axs[1].plot(mass_bins, final_results[:, 2, 1], color=cmap(norm(z)), lw=1, alpha=0.5, linestyle=ls) 
        axs[1].errorbar(mass_bins, final_results[:, 2, 1],
                        yerr=[final_results[:, 2, 1]-final_results[:, 2, 0], 
                              final_results[:, 2, 2]-final_results[:, 2, 1]],
                        color=cmap(norm(z)), fmt='.')
                
simus = ['Hydro', 
         # 'DM'
         ]
for simu in simus:
        
        # Plot TNG300
        for snap in TNG300_snaps:
                TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
                TNG300_z, TNG300_mass_bins, TNG300_final_results = load_stats(TNG300_dir, snap)
                plot_depth(TNG300_z, TNG300_mass_bins, TNG300_final_results, simu)
                plot_width(TNG300_z, TNG300_mass_bins, TNG300_final_results, simu)
                
        # Plot MTNG
        for snap in MTNG_snaps:
                MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
                MTNG_z, MTNG_mass_bins, MTNG_final_results = load_stats(MTNG_dir, snap)
                plot_depth(MTNG_z, MTNG_mass_bins, MTNG_final_results, simu)
                plot_width(MTNG_z, MTNG_mass_bins, MTNG_final_results, simu)

axs[0].set_xscale('log')
axs[0].set_ylabel("Depth")

axs[1].set_xscale('log')
axs[1].set_xlabel('Mass [$M_\\odot$]')
axs[1].set_ylabel("Width")

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig3')
plt.close()