import os
import numpy as np
from func import *
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

print('')
print(f'>>> Plot depth and width vs z (accret) <<<')
print('')

def accret_depth(inputs, a, b, c, e, f):
    zval, Gamma = inputs
    
    return a*Gamma + b*Gamma**2 + c*Gamma**3 + e/(zval+1) + f/Gamma

def accret_width(inputs, a, b, c, d, e, f,):
    zval, Gamma = inputs
    return  a*Gamma + b*Gamma**2 + c*Gamma**3 + d/(zval+1) + f*zval + e

root_dir = 'result/bootstrap_stats/with_accret/'
simus = ['Hydro', 'DM']
features = ['width_dimless' , 
          #  'abs_depth', 
            'depth']

for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = accret_depth
        elif feature == 'width_dimless':
            fit_func = accret_width

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Load data 
        # ============================================================================================

        TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        min_bin, max_bin, bin_width = 1, 6, 1
        num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

        cmap = plt.get_cmap('vanimo', num_bins)
        bound = np.linspace(min_bin, max_bin+0.1, num_bins) 
        norm = BoundaryNorm(bound, cmap.N)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                        ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        cb.set_label(r'$\Gamma$')
        # Reduce colormap ticks sf
        from matplotlib.ticker import FuncFormatter
        def custom_format(x, pos):
            return f'{x:.1f}'  
        cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 

        # ============================================================================================
        # Plot
        # ============================================================================================
        
        all_bin, all_x_med = [], []
        all_y_med, all_y_min, all_y_max = [], [], []
        
        for bin_val in all_bins:
            
            # Plot TNG300
            TNG300_z, TNG300_feat = load_stats_per_bin(TNG300_dir, 'snap_{}_Rsp_stats.npy', TNG300_snaps, 
                                                       'z', feature,
                                                       'accret_bins', bin_val)
            plot_feature_vs_z(simu, bin_val, TNG300_z, TNG300_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(TNG300_z)
            all_y_med.append(TNG300_feat['median'])
            all_y_min.append(TNG300_feat['min'])
            all_y_max.append(TNG300_feat['max'])
            all_bin.append([bin_val]*len(TNG300_z))
                
            # # Plot MTNG
            # MTNG_z, MTNG_feat = load_stats_per_bin(MTNG_dir, 'snap_{}_Rsp_stats.npy', MTNG_snaps,
            #                                        'z', feature,
            #                                        'accret_bins', bin_val)
            # plot_feature_vs_z(simu, bin_val, MTNG_z, MTNG_feat, [axs, cmap, norm])
            
            # all_x_med.append(MTNG_z)    
            # all_y_med.append(MTNG_feat['median'])
            # all_y_min.append(MTNG_feat['min'])
            # all_y_max.append(MTNG_feat['max'])
            # all_bin.append([bin_val]*len(MTNG_z))
            
        # Final edit
        axs.set_xlabel('z')
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        
        # ============================================================================================
        # Fitting
        # ============================================================================================
        
        # Concatenate the data to one dimension
        all_x_med = np.concatenate(all_x_med)
            
        all_y_med = np.concatenate(all_y_med)
        all_y_min = np.concatenate(all_y_min)
        all_y_max = np.concatenate(all_y_max)

        all_bin = np.concatenate(all_bin)
        
        # Fit the data
        popt, red_chi2, y_fit, axs = fitting(all_bin, 
                                             all_x_med,
                                              all_y_med, all_y_min, all_y_max,
                                              fit_func, [axs, cmap, norm])
        print(f'{simu} {feature} popt: {popt}')
        print(f'{simu} {feature} red_chi2: {red_chi2}')
        axs.legend()
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig4_accret_{simu}_{feature}')
        plt.close()