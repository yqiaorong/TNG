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

X_d, Xerr_d, Y_d, Yerr_d, z_d = [], [], [], [], []
x, xerr, z, y, yerr = plot_data(TNG300_dir, TNG300_snaps, [TNG300_z, TNG300_mass_cuts, TNG300_accret_med, TNG300_accret_width], 
                                axs, cmap, norm, feat='depth')
X_d.append(x)
Xerr_d.append(xerr)
Y_d.append(y)
Yerr_d.append(yerr)
z_d.append(z)
x, xerr, z, y, yerr = plot_data(MTNG_dir, MTNG_snaps, [MTNG_z, MTNG_mass_cuts, MTNG_accret_med, MTNG_accret_width], 
                                axs, cmap, norm, feat='depth')
X_d.append(x)
Xerr_d.append(xerr)
Y_d.append(y)
Yerr_d.append(yerr)
z_d.append(z)


X_d = np.concatenate(X_d)
Xerr_d = np.concatenate(Xerr_d)
Y_d = np.concatenate(Y_d)
Yerr_d = np.concatenate(Yerr_d, axis=1)

# # add 0.01 to every element in Yerr_d
# Yerr_d = [Yerr_d[i] + 0 for i in range(len(Yerr_d))]

z_d = np.concatenate(z_d)
# print(X_d.shape, Xerr_d.shape, Y_d.shape, Yerr_d.shape, z_d.shape)

save_dict = {'accret_rate': X_d,
             'accret_rate_err': Xerr_d, 
             'depth': Y_d,
             'depth_err': Yerr_d,
             'z': z_d}
np.save('depth_data.npy', save_dict)

# ============================================================================================
# Fitting 
# ============================================================================================

# # if the array in x shares the same z value, then concatenate the array
# DoF = 1
# poly_fit(X_d, Y_d, z_d, axs, cmap, norm, DoF)

# Fit two params, y(x, z) at the same time
x_z_fit_d(X_d, z_d, Y_d, Yerr_d, axs, cmap, norm)

# Final edit
axs.set_ylabel(r"$\mathcal{D}$")
axs.set_xlabel(r'$\Gamma$')
axs.legend()

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig6_fit_D.png')
plt.close()