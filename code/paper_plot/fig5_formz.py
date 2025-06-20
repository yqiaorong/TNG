import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import BoundaryNorm
from func import load_stats, plot_feature, fitting

# ============================ CONFIG ============================
simu = 'Hydro'
x_type = 'formzOLD'
features = ['depth', 'width_dimless']
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$"]
xlabel = r'$z_{\rm form}$'
save_path = f'result/paper_plots/'

# Input directories
fig4_root = f'result/bootstrap_stats_DK14/with_{x_type}_perMassCut/'  # Right column (mass)
fig3_root = f'result/bootstrap_stats_DK14/with_{x_type}/'              # Left column (z)

# Prepare mass bins and color map
min_bin, max_bin, bin_width = 13, 15.5, 0.5
num_bins = int((max_bin - min_bin)/bin_width)

cmap_mass = plt.get_cmap('plasma', num_bins)
all_bins = np.linspace(min_bin, max_bin, num_bins+1)
bound_mass = np.logspace(min_bin, max_bin, num_bins+1)
norm_mass = BoundaryNorm(bound_mass, cmap_mass.N)

cmap_z = plt.get_cmap('viridis')
all_z = np.array([0.0, 1.0, 2.0])
bound_z = np.array([0.0, 0.5, 1.5, 2.5])
norm_z = BoundaryNorm(all_z, cmap_z.N)

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================

def depth(inputs, a, b):
    x, zval = inputs
    return a*np.log10(x) + b

def depth_from_mass(inputs, A, a, B):
    x, z = inputs
    return A * (a-np.log10(x))**B

def width(inputs, a, b):
    x, zval = inputs
    return a*np.log(x) + b

def width_from_mass(inputs, A, a, B):
    x, z = inputs
    return A * (a-np.log10(x))**B

# ============================ SETUP FIGURE ============================
# fig, axs = plt.subplots(2, 2, figsize=(9, 8), dpi=500, sharex=True, sharey='row', constrained_layout=True,
#                         # gridspec_kw={'height_ratios': [1, 1], 'width_ratios': [1, 1]}
#                         )

import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(7, 6), dpi=500)
gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 0.05], hspace=0.2)  # 2 rows for plots, 1 for colorbar

axs = np.empty((2, 2), dtype=object)  # 2 rows × 3 columns subplot array

# # # Create subplots for top 2 rows
# for col in range(2):
#     # axs[row, col] = fig.add_subplot(gs[row, col])
#     axs[0, col] = fig.add_subplot(gs[0, col], )
#     axs[1, col] = fig.add_subplot(gs[1, col], sharex=axs[0, col]) 
#     axs[0, col].tick_params(labelbottom=False)
    
ax10 = fig.add_subplot(gs[1, 0])
ax00 = fig.add_subplot(gs[0, 0], sharex=ax10)

ax11 = fig.add_subplot(gs[1, 1], sharey=ax10)
ax01 = fig.add_subplot(gs[0, 1], sharex=ax11, sharey=ax00)

axs = np.array([[ax00, ax01],
                [ax10, ax11]])

    
# ax00.xaxis.set_visible(False)
ax01.tick_params(labelleft=False)
ax11.tick_params(labelleft=False)
ax00.tick_params(labelbottom=False)
ax01.tick_params(labelbottom=False)

# ============================ PLOT EACH FEATURE ============================
for ifeat, feature in enumerate(features):
    # ---------- Left Column: snapshot z (fig3 style) ----------
    all_z_vals, x_med, x_min, x_max = [], [], [], []
    y_med, y_min, y_max = [], [], []
    for snap in [264]:
        zval, xdata, ydata = load_stats(os.path.join(fig3_root, 'MTNG', f'{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'),
                                        f'snap_{snap}_Rsp_stats.npy', f'med_{x_type}', feature)
        mask = xdata['median'] >= 0.25
        xdata = {k: v[mask] for k, v in xdata.items()}
        ydata = {k: v[mask] for k, v in ydata.items()}
        plot_feature(simu, np.round(zval, 1), xdata, ydata, [axs[ifeat, 0], cmap_z, norm_z])
        x_med.append(xdata['median'])
        x_min.append(xdata['min'])
        x_max.append(xdata['max'])
        y_med.append(ydata['median'])
        y_min.append(ydata['min'])
        y_max.append(ydata['max'])
        all_z_vals.append([zval]*len(xdata['median']))

    x_med = np.concatenate(x_med)
    x_min = np.concatenate(x_min)
    x_max = np.concatenate(x_max)
    y_med = np.concatenate(y_med)
    y_min = np.concatenate(y_min)
    y_max = np.concatenate(y_max)
    all_z_vals = np.concatenate(all_z_vals)

    valid = ~np.isnan(y_med)
    x_med, all_z_vals, y_med, y_min, y_max = x_med[valid], all_z_vals[valid], y_med[valid], y_min[valid], y_max[valid]

    popt, _, _, _, _, _ = fitting(depth if feature=='depth' else width, 
                                  [x_med, all_z_vals, y_med, y_min, y_max],
                                  [x_type, 'z', feature], [axs[ifeat, 0], cmap_z, norm_z, '--'], bootstrap=True)

    popt, _, _, _, _, _ = fitting(depth_from_mass if feature=='depth' else width_from_mass, 
                                  [x_med, all_z_vals, y_med, y_min, y_max],
                                  [x_type, 'z', feature], [axs[ifeat, 0], cmap_z, norm_z, 'dotted'], bootstrap=False)

    axs[ifeat, 0].set_ylabel(Ylabels[ifeat])

    # ---------- Right Column: perMassCut (fig5 style) ----------
    all_x2, x_med, x_min, x_max = [], [], [], []
    y_med, y_min, y_max = [], [], []
    mass_dir = os.path.join(fig4_root, 'MTNG', f'{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/')
    for fname in os.listdir(mass_dir):
        zval, xdata, ydata = load_stats(mass_dir, fname, f'med_{x_type}', feature)
        bin_val = 10**(10 + float(fname.split('_')[3])/10)
        mask = xdata['median'] >= 0.25
        xdata = {k: v[mask] for k, v in xdata.items()}
        ydata = {k: v[mask] for k, v in ydata.items()}
        plot_feature(simu, bin_val, xdata, ydata, [axs[ifeat, 1], cmap_mass, norm_mass])
        x_med.append(xdata['median'])
        x_min.append(xdata['min'])
        x_max.append(xdata['max'])
        y_med.append(ydata['median'])
        y_min.append(ydata['min'])
        y_max.append(ydata['max'])
        all_x2.append(np.repeat(bin_val, len(xdata['median'])))

    x_med = np.concatenate(x_med)
    x_min = np.concatenate(x_min)
    x_max = np.concatenate(x_max)
    y_med = np.concatenate(y_med)
    y_min = np.concatenate(y_min)
    y_max = np.concatenate(y_max)
    all_x2 = np.concatenate(all_x2)

    valid = ~np.isnan(y_med)
    x_med, all_x2, y_med, y_min, y_max = x_med[valid], all_x2[valid], y_med[valid], y_min[valid], y_max[valid]

    popt, _, _, _, _, _ = fitting(depth if feature=='depth' else width, 
                                  [x_med, all_x2, y_med, y_min, y_max],
                                  [x_type, 'mass', feature], [axs[ifeat, 1], None, norm_mass, '--'], bootstrap=True)

    popt, _, _, _, _, _ = fitting(depth_from_mass if feature=='depth' else width_from_mass, 
                                  [x_med, all_x2, y_med, y_min, y_max],
                                  [x_type, 'mass', feature], [axs[ifeat, 1], None, norm_mass, 'dotted'], bootstrap=False)

# ============================ LABELING & COLORBARS ============================
axs[1, 0].set_xlabel(xlabel)
axs[1, 1].set_xlabel(xlabel)

# Create long horizontal colorbar in the third row
cbar_ax = fig.add_subplot(gs[2, 1]) 
cb = fig.colorbar(cm.ScalarMappable(norm=norm_mass, cmap=cmap_mass), cax=cbar_ax,
                  orientation='horizontal', spacing='proportional', ticks=bound_mass)
cb.set_label(r'$M_{200m}/M_{\odot}$')
cb._set_scale('log')

fig.subplots_adjust(
    left=0.1, right=0.98,
    top=0.97,  bottom=0.05,   
    wspace=0.02, 
)

shift = 0.03
for ax in (ax10, ax11):
    pos = ax.get_position()     
    print(pos)# [x0, y0, width, height]
    ax.set_position([pos.x0, pos.y0 + shift,        # y0 往上挪
                     pos.width, pos.height])
    print(ax.get_position())
    print('')

fig.canvas.draw()

# ============================ SAVE ============================
os.makedirs(save_path, exist_ok=True)
# fig.subplots_adjust(bottom=0.15) 
plt.savefig(save_path+f'fig5_{x_type}_{simu}')
plt.close()