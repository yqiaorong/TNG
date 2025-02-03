""""This script plots width as a function of accretion rate.
    Hydro simulation only."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm


plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot width vs accretion rate <<<')
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

fig, axs = plt.subplots(1, 1, figsize = (5, 4), dpi=500, sharex=True, constrained_layout=True)
cmap = plt.get_cmap('managua', num_z)
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
# Load accretion rate width
# ============================================================================================

_, _, _, TNG300_accret_width = load_accret(f'result/accretion_rate_plot/TNG300/sim_205_1250_Hydro/TNG300_Hydro_accret_stats.npy', 
                                                                 width='percentile')

_, _, _, MTNG_accret_width = load_accret(f'result/accretion_rate_plot/MTNG/Hydro-Arepo/MTNG_Hydro_accret_stats.npy', 
                                                           width='percentile')

# ============================================================================================
# Plot 
# ============================================================================================

X_w, Xerr_w, Y_w, Yerr_w, z_w = [], [], [], [], []
x, xerr, z, y, yerr = plot_data(TNG300_dir, TNG300_snaps, [TNG300_z, TNG300_mass_cuts, TNG300_accret_med, TNG300_accret_width], 
                                axs, cmap, norm, feat='width')
X_w.append(x)
Xerr_w.append(xerr)
Y_w.append(y)
Yerr_w.append(yerr)
z_w.append(z)
x, xerr, z, y, yerr = plot_data(MTNG_dir, MTNG_snaps, [MTNG_z, MTNG_mass_cuts, MTNG_accret_med, MTNG_accret_width], 
                                axs, cmap, norm, feat='width')
X_w.append(x)
Xerr_w.append(xerr)
Y_w.append(y)
Yerr_w.append(yerr)
z_w.append(z)

X_w = np.concatenate(X_w)
Xerr_w = np.concatenate(Xerr_w)
Y_w = np.concatenate(Y_w)
Yerr_w = np.concatenate(Yerr_w, axis=1)

z_w = np.concatenate(z_w)

# ============================================================================================
# Fitting 
# ============================================================================================

# if the array in x shares the same z value, then concatenate the array
# poly_fit(X_w, Y_w, Yerr_w, z_w, axs, cmap, norm, 4)

# Fit two params, y(x, z) at the same time
x_z_fit_w(X_w, z_w, Y_w, Yerr_w, axs, cmap, norm)
# x_z_fit_w(X_w[z_w == 0], z, Y_w[z_w == 0], Yerr_w[:, z_w == 0], axs)

# Final edit
axs.set_ylabel(r"$\mathcal{W}$")
axs.set_xlabel(r'$\Gamma$')
axs.legend()

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig6_fit_W.png')
plt.close()