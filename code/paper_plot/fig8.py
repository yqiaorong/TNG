""""This script plots the absolute depth as a function of mass for a few redshifts.
    Both Hydro and DM-only simulations."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot absolute depth vs mass for a few redshifts <<<')
print('')

root_dir = 'result/bootstrap_stats_phys/'

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

# ============================================================================================
# Set up the colorbar
# ============================================================================================

TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

all_z = load_all_z(TNG300_dir, TNG300_snaps)

MTNG_snaps = [264, 237, 214, 179, 151, 129]
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)

# Set up the colorbar
# -----------------------------------------------------------------------------------------
# from matplotlib.colors import LinearSegmentedColormap
# red_list = [# '#EE9D9F', 
#             '#DE6A69', '#C84747', '#982B2D','#6A0624', '#3D011A'] # light to dark
# blue_list = ['#89CAEA','#4596CD', '#0B75B3', '#015696', '#012A61']
# cmap = LinearSegmentedColormap.from_list('my_cmap', red_list)
cmap = plt.get_cmap('managua', len(TNG300_snaps))
bound = all_z
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# -----------------------------------------------------------------------------------------

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.ax.tick_params(axis='x', rotation=70) 
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

def plot_abs_depth(z, mass_bins, final_results, simu):
        if simu == 'DM':
                ls = '--'
        else:
                ls = '-'
        
        # -----------------------------------------------------------------------------------------
        axs.plot(mass_bins, final_results[:, 4, 1], color=cmap(norm(z)), lw=1, alpha=0.5, linestyle=ls) 
        axs.errorbar(mass_bins, final_results[:, 4, 1],
                        yerr=[final_results[:, 4, 1]-final_results[:, 4, 0], 
                              final_results[:, 4, 2]-final_results[:, 4, 1]],
                        color=cmap(norm(z)), fmt='.')
        # axs.plot(mass_bins, abs(final_results[:, 4, 1]), color=cmap(norm(z)), lw=1, alpha=0.5, linestyle=ls) 
        # axs.errorbar(mass_bins, abs(final_results[:, 4, 1]),
        #                 yerr=[abs(final_results[:, 4, 1]-final_results[:, 4, 0]), 
        #                       abs(final_results[:, 4, 2]-final_results[:, 4, 1])],
        #                 color=cmap(norm(z)), fmt='.')
        # -----------------------------------------------------------------------------------------
                
simus = ['Hydro', 
         # 'DM'
         ]
for simu in simus:
        
        # Plot TNG300
        for snap in TNG300_snaps:
                TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
                TNG300_z, TNG300_mass_bins, TNG300_final_results = load_stats(TNG300_dir, snap)
                plot_abs_depth(TNG300_z, TNG300_mass_bins, TNG300_final_results, simu)
                
        # Plot MTNG
        for snap in MTNG_snaps:
                MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
                MTNG_z, MTNG_mass_bins, MTNG_final_results = load_stats(MTNG_dir, snap)
                plot_abs_depth(MTNG_z, MTNG_mass_bins, MTNG_final_results, simu)

axs.set_xscale('log')
axs.set_ylabel(r"abs $\mathcal{D}$")
axs.set_xlabel('Mass [$M_\\odot$]')

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig8')
plt.close()