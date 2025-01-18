""""This script plots depth as a function of accretion rate.
    Hydro simulation only."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm


plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot depth vs accretion rate <<<')
print('')

root_dir = 'result/bootstrap_stats_phys/'

# ============================================================================================
# Load z
# ============================================================================================

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_snaps = [99, 67, 40, 25, 13, 8]
TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

all_z = load_all_z(TNG300_dir, TNG300_snaps)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_snaps = [264, 214, 151]
MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)
num_z = len(TNG300_snaps)

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(2, 1, figsize = (4, 6), dpi=500, sharex=True, constrained_layout=True)

from matplotlib.colors import LinearSegmentedColormap
# red_list = [# '#EE9D9F', '#DE6A69', 
#             '#C84747', '#982B2D','#6A0624', '#3D011A'] # light to dark
# blue_list = ['#89CAEA','#4596CD', '#0B75B3', '#015696', '#012A61']
# cmap = LinearSegmentedColormap.from_list('my_cmap', red_list)
cmap = plt.get_cmap('vanimo', num_z)
bound = all_z
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.set_label('z')

# ============================================================================================
# Load accretion rate
# ============================================================================================

TNG300_mass_cuts, TNG300_z, TNG300_accret_med = load_accret(f'result/accretion_rate_plot/TNG300/sim_205_1250_Hydro/TNG300_Hydro_accret_stats.npy')

MTNG_mass_cuts, MTNG_z, MTNG_accret_med = load_accret(f'result/accretion_rate_plot/MTNG/Hydro-Arepo/MTNG_Hydro_accret_stats.npy')

# ============================================================================================
# Plot 
# ============================================================================================

_, _ = plot_data(TNG300_dir, TNG300_snaps, [TNG300_z, TNG300_mass_cuts, TNG300_accret_med], axs[0], cmap, norm, feat='depth')
_, _ = plot_data(MTNG_dir, MTNG_snaps, [MTNG_z, MTNG_mass_cuts, MTNG_accret_med], axs[0], cmap, norm, feat='depth')

_, _ = plot_data(TNG300_dir, TNG300_snaps, [TNG300_z, TNG300_mass_cuts, TNG300_accret_med], axs[1], cmap, norm, feat='width')
_, _ = plot_data(MTNG_dir, MTNG_snaps, [MTNG_z, MTNG_mass_cuts, MTNG_accret_med], axs[1], cmap, norm, feat='width')

# Final edit
axs[0].set_ylabel(r"$\mathcal{D}$")
axs[1].set_ylabel(r"$\mathcal{W}$")
axs[1].set_xlabel(r'$\Gamma$')

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig6.png')
plt.close()