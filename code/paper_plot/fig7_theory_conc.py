import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *


bin_type = 'NFWconc'
xlabel = r'$c$'
theory_types = ['rho_s', 'r_s', 'alpha', 'r_t', 'beta', 'gamma', 'rho_g', 'b_e', 'S_e', 'R200']
    
print('')
print(f'>>> Plot depth and width vs {bin_type} _perMassCut <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}_perMassCut/'
# root_dir2 = f'result/bootstrap_stats2/with_{bin_type}_perMassCut/'
simus = ['Hydro']
features = ['abs_depth', 'depth', 'width_dimless', 'DWratio']

save_dir = f'result/paper_plots/fig7/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


for simu in simus:
    for feature in features:
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

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
                          ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        
        # Reduce colormap ticks sf
        from matplotlib.ticker import FuncFormatter
        def custom_format(x, pos):
            return f'{x:.1f}'  
        cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        cb.ax.tick_params(axis='x', rotation=0) 
        
        cb.set_label(r'$M_{\odot}$')
        cb._set_scale('log')
            
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
            plot_theory_feature(bin_val, MTNG_bin_data_acc['median'][sort_idx], np.array(theory_feats)[sort_idx], 
                                [axs, cmap, norm], label=f'z=0.0')
                    
            plot_feature(simu, bin_val, MTNG_bin_data_acc, MTNG_feat_acc, [axs, cmap, norm], label=f'z={MTNG_z:.1f}')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        if feature == 'DWratio':
            axs.set_ylim(0.5, 4)
        if bin_type == 'mass':
            axs.set_xscale('log')
        axs.set_xlabel(xlabel)
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'abs_depth':
            Y_label = r"-|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        # axs.legend(loc='best')

        plt.savefig(f'{save_dir}/{bin_type}_{simu}_{feature}')
        plt.close()