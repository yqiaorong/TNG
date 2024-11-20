""""This script plots depth (top panel) and width (bottom panel) as a function of mass at z = 0.
    Hydro-simulations only."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

print('')
print('>>> Plot Depth and width vs mass at z = 0 <<<')
print('')

# ============================================================================================
# Load data
# ============================================================================================

root_dir = 'result/bootstrap_stats_phys/'

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_snaps = 99
TNG300_data = np.load(f'{TNG300_dir}/snap_{TNG300_snaps}_Rsp_stats.npy', allow_pickle=True).item()

TNG300_z             = np.round(TNG300_data['z'], 3)
TNG300_mass_bins     = TNG300_data['mass_bins']
TNG300_final_results = TNG300_data['final_results']

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_snaps = 264
MTNG_data = np.load(f'{MTNG_dir}//snap_{MTNG_snaps}_Rsp_stats.npy', allow_pickle=True).item()

MTNG_z             = np.round(MTNG_data['z'], 3)
MTNG_mass_bins     = MTNG_data['mass_bins']
MTNG_final_results = MTNG_data['final_results']

print(TNG300_z, MTNG_z)

min_mass = min(min(TNG300_mass_bins), min(MTNG_mass_bins))
max_mass = max(max(TNG300_mass_bins), max(MTNG_mass_bins))
print(min_mass, max_mass)
num_mass_bins = int((max_mass - min_mass) / 0.5)
print(num_mass_bins)

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(2, 1, figsize=(4, 6), dpi=500, sharex=True)

from matplotlib.colors import LinearSegmentedColormap
red_list = [# '#EE9D9F', '#DE6A69', 
            '#C84747', '#982B2D','#6A0624', '#3D011A'] 
blue_list = ['#89CAEA','#4596CD','#0B75B3','#015696','#012A61','#053061', ]# light to dark
cmap = LinearSegmentedColormap.from_list('my_cmap', blue_list)
# cmap = plt.get_cmap('plasma', num_mass_bins)

bound = [10**(10+mass) for mass in np.linspace(min_mass, max_mass, num_mass_bins)]
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label('Mass [$M_\\odot$]')
cb._set_scale('log')

# ============================================================================================
# Plot Depth
# ============================================================================================

# Plot TNG300
for imass, mass in enumerate(TNG300_mass_bins):
        
        axs[0].plot(10**(10+mass), TNG300_final_results[imass, 1, 1], color=cmap(norm(10**(10+mass))))
        if imass < len(TNG300_mass_bins)-1:
                axs[0].plot([10**(10+TNG300_mass_bins[imass]), 10**(10+TNG300_mass_bins[imass+1])], 
                            [TNG300_final_results[imass, 1, 1], TNG300_final_results[imass+1, 1, 1]], 
                             color=cmap(norm(10**(10+TNG300_mass_bins[imass]))), lw=1) 
        axs[0].errorbar(10**(10+mass), TNG300_final_results[imass, 1, 1],
                        yerr=[[TNG300_final_results[imass, 1, 1]-TNG300_final_results[imass, 1, 0]], 
                                [TNG300_final_results[imass, 1, 2]-TNG300_final_results[imass, 1, 1]]],
                        color=cmap(norm(10**(10+mass))), fmt='.')

# Plot MTNG
for imass, mass in enumerate(MTNG_mass_bins):
        
        axs[0].plot(10**(10+mass), MTNG_final_results[imass, 1, 1], color=cmap(norm(10**(10+mass))))
        if imass < len(MTNG_mass_bins)-1:
                axs[0].plot([10**(10+MTNG_mass_bins[imass]), 10**(10+MTNG_mass_bins[imass+1])], 
                            [MTNG_final_results[imass, 1, 1], MTNG_final_results[imass+1, 1, 1]],
                             color=cmap(norm(10**(10+MTNG_mass_bins[imass]))), lw=1)
        axs[0].errorbar(10**(10+mass), MTNG_final_results[imass, 1, 1],
                        yerr=[[MTNG_final_results[imass, 1, 1]-MTNG_final_results[imass, 1, 0]],
                                [MTNG_final_results[imass, 1, 2]-MTNG_final_results[imass, 1, 1]]],
                        color=cmap(norm(10**(10+mass))), fmt='.')

axs[0].set_xscale('log')
axs[0].set_ylabel("Depth")
axs[0].set_xlim(10**(10+min_mass), 10**(10+max_mass))

# ============================================================================================
# Plot Width
# ============================================================================================

# Plot TNG300
for imass, mass in enumerate(TNG300_mass_bins):
                
        axs[1].plot(10**(10+mass), TNG300_final_results[imass, 2, 1], color=cmap(norm(10**(10+mass))))
        if imass < len(TNG300_mass_bins)-1:
                axs[1].plot([10**(10+TNG300_mass_bins[imass]), 10**(10+TNG300_mass_bins[imass+1])], 
                            [TNG300_final_results[imass, 2, 1], TNG300_final_results[imass+1, 2, 1]],
                        color=cmap(norm(10**(10+TNG300_mass_bins[imass]))), lw=1)
        axs[1].errorbar(10**(10+mass), TNG300_final_results[imass, 2, 1],
                        yerr=[[TNG300_final_results[imass, 2, 1]-TNG300_final_results[imass, 2, 0]],
                                [TNG300_final_results[imass, 2, 2]-TNG300_final_results[imass, 2, 1]]],
                        color=cmap(norm(10**(10+mass))), fmt='.')
        
# Plot MTNG
for imass, mass in enumerate(MTNG_mass_bins):
                        
        axs[1].plot(10**(10+mass), MTNG_final_results[imass, 2, 1], color=cmap(norm(10**(10+mass))))
        if imass < len(MTNG_mass_bins)-1:
                axs[1].plot([10**(10+MTNG_mass_bins[imass]), 10**(10+MTNG_mass_bins[imass+1])], 
                            [MTNG_final_results[imass, 2, 1], MTNG_final_results[imass+1, 2, 1]],
                        color=cmap(norm(10**(10+MTNG_mass_bins[imass]))), lw=1)
        axs[1].errorbar(10**(10+mass), MTNG_final_results[imass, 2, 1],
                        yerr=[[MTNG_final_results[imass, 2, 1]-MTNG_final_results[imass, 2, 0]],
                                [MTNG_final_results[imass, 2, 2]-MTNG_final_results[imass, 2, 1]]],
                        color=cmap(norm(10**(10+mass))), fmt='.')

axs[1].set_xscale('log')
axs[1].set_xlabel('Mass [$M_\\odot$]')
axs[1].set_ylabel("Width")
axs[1].set_xlim(10**(10+min_mass), 10**(10+max_mass))

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig2.png')
plt.close()