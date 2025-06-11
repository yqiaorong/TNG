import os
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib import pyplot as plt
plt.style.use('code/style.mplstyle')


sim = 'MTNG/Hydro-Arepo/MTNG-L500-4320-A/'
snapnum, z = 129, 2.0
xval_name = 'accretions' 
yval_name = 'mass'

if xval_name == 'formzOLD':
    xlabel = r'$z_{\rm form}$'
elif xval_name == 'formzSub':
    xlabel = r'$z_{\rm form}$ (half subhalo mass)'
elif xval_name == 'accretions':
    xlabel = r'$\Gamma$'
elif xval_name == 'peakHeight':
    xlabel=r'$v$'
    
if yval_name == 'NFWconc':
    ylabel = r'$c$'
elif yval_name == 'mass':
    ylabel = r'$M_{200m} / M_{\odot}$'
    

# Plot the accretion rate vs mass
fig, axs = plt.subplots(1, 1, figsize=(4, 3.3), dpi=500)

# Set up the colorbar
num_bins = 10
x_start, x_end = 0, 2.5
x_width = (x_end - x_start) / (num_bins-1)
xval_bins = np.linspace(x_start, x_end, num_bins)

cmap = plt.get_cmap('plasma', num_bins)
norm = BoundaryNorm(xval_bins, cmap.N)

halos_dir = f'result/DMhalo_density_profiles/{sim}/snap_{snapnum}/final_densities/'
halos_fnames = os.listdir(halos_dir)
print(halos_fnames)
for fname in halos_fnames:
    data = np.load(halos_dir+fname, allow_pickle=True).item()
    print(data.keys())
    
    if yval_name == 'mass':
        yvals = data['halo_M_Mean200']
        yvals = 10**10 * yvals  # convert to 10^10 Msun
    else:
        yvals = data[yval_name]         
    xvals = data[xval_name]    
    
    # Select halos 
    stats_yvals = np.zeros((len(xval_bins)-1, 3))
    for idx, i in enumerate(xval_bins[:-1]):
        mask = (xvals >= i) & (xvals < i+x_width)
        selected_yvals = yvals[mask]
        selected_vals = xvals[mask]
        # Plot
        axs.scatter(selected_vals, selected_yvals,  s=1, color=cmap(norm(i)), alpha=0.5)
        # Calculate the median and 16th/84th percentiles
        stats_yvals[idx, 0] = np.median(selected_yvals)
        stats_yvals[idx, 1] = np.percentile(selected_yvals, 16)
        stats_yvals[idx, 2] = np.percentile(selected_yvals, 84)
    
    # Plot the median and percentiles
    axs.errorbar(xval_bins[:-1] + x_width/2, stats_yvals[:, 0], 
                yerr=[stats_yvals[:, 0] - stats_yvals[:, 1], stats_yvals[:, 2] - stats_yvals[:, 0]], 
                fmt='.', linestyle='-', color='black')
    
    # Place text
    axs.text(0.85, 0.9, f'z={z}', transform=axs.transAxes, fontsize=12)
if yval_name == 'mass':
    axs.set_ylim(10**13, 10**15.5)
    axs.set_yscale('log')
elif yval_name == 'NFWconc':
    axs.set_ylim(0, 50)
axs.set_xlabel(xlabel)
axs.set_ylabel(ylabel)
# plt.tight_layout()
save_dir = f'result/paper_plots/fig_vs_{yval_name}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/{xval_name}_{snapnum}')
plt.close()