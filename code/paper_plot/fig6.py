import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, get_bins, plot_feature

bin_type = 'mergerz'
if bin_type == 'peakHeight':
    xlabel=r'$v$'
elif bin_type == 'NFWconc':
    xlabel = r'$c$'
elif bin_type == 'mergerz':
    xlabel = r'$z_{\rm merger}$'
elif bin_type == 'formz':
    xlabel = r'$z_{\rm form}$'
elif bin_type == 'formzOLD':
    xlabel = r'$z_{\rm form}$ (half mass)'
elif bin_type == 'formzSub':
    xlabel = r'$z_{\rm form}$ (half subhalo mass)'
# elif bin_type in ['accretions', 'accretionsOLD', 'accret']:
#     xlabel = r'$\Gamma$'
# elif bin_type == 'mass':
#     xlabel = r'$M_{\odot}$'
    
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}_perMassCut/'
simus = ['Hydro']
features = ['width_dimless','depth', 'DWratio']


for simu in simus:
    for feature in features:
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.3), dpi=500, sharex=True, constrained_layout=True)

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
        # Plot
        # ============================================================================================
              
        # # Plot TNG300       
        # for snap in TNG300_snaps:
        #     TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
        #                                                     f'med_{bin_type}', feature)
        #     plot_feature(simu, TNG300_z, TNG300_bin_data, TNG300_feat, [axs, cmap, norm])
            
        # Plot MTNG
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        MTNG_list = os.listdir(MTNG_dir)
        for fname in MTNG_list:
            print(fname)
            _, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, fname, f'med_{bin_type}', feature)    
            bin_val = 10**(10+float(fname.split('_')[3])/10)
            plot_feature(simu, bin_val, MTNG_bin_data, MTNG_feat, [axs, cmap, norm])

        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xlabel(xlabel)
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)

        # axs.legend(loc='best')
        save_dir = f'result/paper_plots/fig6/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig6_{bin_type}_perMassCut_{simu}_{feature}')
        plt.close()