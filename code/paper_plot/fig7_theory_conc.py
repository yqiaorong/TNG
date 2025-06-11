import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *
from matplotlib.ticker import FormatStrFormatter


bin_type = 'NFWconc'
xlabel = r'$c$'
theory_types = ['rho_s', 'r_s', 'alpha', 'r_t', 'beta', 'gamma', 'rho_g', 'b_e', 'S_e', 'R200']
    
print('')
print(f'>>> Plot depth and width vs {bin_type} _perMassCut <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}_perMassCut/'
simu = 'Hydro'
features = ['depth', 'width_dimless', 'DWratio']
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$", r"$\mathcal{D}/\mathcal{W}$"]



# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(3, 1, figsize=(4, 9), dpi=500, sharex=True, constrained_layout=True)

# ============================================================================================
# Set up the colorbar
# ============================================================================================

# Use mass cuts
min_bin, max_bin, bin_width = 13, 15, 0.25
num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

cmap = plt.get_cmap('plasma', num_bins)
bound = np.logspace(min_bin, max_bin+0.01, num_bins+1) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=axs[-1], orientation='horizontal', spacing='proportional', ticks=bound)

# Reduce colormap ticks sf
# from matplotlib.ticker import FuncFormatter
# def custom_format(x, pos):
#     return f'{x:.1f}'  
# cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
# cb.ax.tick_params(axis='x', rotation=0) 

cb.set_label(r'$M_{200m}/M_{\odot}$')
cb._set_scale('log')

save_dir = f'result/paper_plots/fig7/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)



for ifeat, feature in enumerate(features):
    axs[ifeat].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    
    # ============================================================================================
    # Plot the theory
    # ============================================================================================
    
    radius = np.logspace(-2, np.log10(5), 1024)
    
    # Plot MTNG
    MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
    MTNG_list = os.listdir(MTNG_dir)
    for fname in MTNG_list:
        print(fname)
        
        # Load conc data
        MTNG_z, MTNG_bin_data_acc, MTNG_feat_acc = load_stats(f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/', 
                                                        fname, f'med_{bin_type}', feature) 
        bin_val = 10**(10+float(fname.split('_')[3])/10)
        MTNG_params = load_stats_key(f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/', fname, 'params') 
        
        # Assemble params
        params = []
        for key, val in MTNG_params.items():
            params.append(val[1])
        params = np.array(params)

        theory_feats = []
        for i in range(params.shape[1]):
            slopes = density_gradient_profile(radius, *params[:,i])
            
            # depth
            min_grad = np.min(slopes)
            min_grad_idx = np.argmin(slopes)
            if min_grad_idx > 0:
                
                plot_slope_profile(radius, slopes, save_dir, f'profile_{fname.split('.')[0][:-10]}_{i}')
                
                print(min_grad, min_grad_idx)
                left_data = slopes[:min_grad_idx]
                right_data = slopes[min_grad_idx:]
                
                max_grad = np.max(right_data)
                depth = max_grad - min_grad
                
                # Width
                half_grad = min_grad + depth/2
                left_idx = np.argmin(np.abs(left_data - half_grad))
                right_idx = min_grad_idx + np.argmin(np.abs(right_data - half_grad))
                width_dimless = radius[right_idx] - radius[left_idx]
                
                # Depth vs width
                DWratio = depth / width_dimless
                
                if feature == 'depth':
                    theory_feats.append(depth)
                elif feature == 'width_dimless':
                    theory_feats.append(width_dimless)
                elif feature == 'DWratio':
                    theory_feats.append(DWratio)
                elif feature == 'abs_depth':
                    theory_feats.append(min_grad)
            else:
                theory_feats.append(np.nan)

        # Sort data
        sort_idx = np.argsort(MTNG_bin_data_acc['median'])
        if feature != 'depth':
            plot_theory_feature(bin_val, MTNG_bin_data_acc['median'][sort_idx], np.array(theory_feats)[sort_idx], 
                            [axs[ifeat], cmap, norm], label=f'z=0.0')
                
        plot_feature(simu, bin_val, MTNG_bin_data_acc, MTNG_feat_acc, 
                     [axs[ifeat], cmap, norm], label=f'z={MTNG_z:.1f}')
        
    if bin_type == 'mass':
        axs[ifeat].set_xscale('log')
    axs[ifeat].set_ylabel(Ylabels[ifeat])
    
# ============================================================================================
# Save the plot
# ============================================================================================
# if feature == 'DWratio':
#     axs.set_ylim(0.5, 4)

axs[2].set_xlabel(xlabel)

plt.savefig(f'{save_dir}/{bin_type}_{simu}')
plt.close()