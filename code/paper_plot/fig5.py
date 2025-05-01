import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

    
print('')
print(f'>>> Plot depth and width vs formzOLD and formz <<<')
print('')


simus = ['Hydro']
features = ['width_dimless','depth']


for simu in simus:
    for feature in features:
        
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)
        
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        all_z = np.linspace(0, 3, 10)
        cmap = plt.get_cmap('viridis', len(all_z))
        if len(all_z) == 1:
            bound = [all_z[0], all_z[0]+0.1]
        else:
            bound = all_z
        norm = BoundaryNorm(bound, cmap.N)
        
        # cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
        #                 ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        # -----------------------------------------------------------------------------------------

        # # Reduce colormap ticks sf
        # from matplotlib.ticker import FuncFormatter
        # def custom_format(x, pos):
        #     return f'{x:.1f}'  
        # cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        # cb.ax.tick_params(axis='x', rotation=70) 
        # cb.set_label('z')
        
        for bin_type, label, ls in zip(['formz', 'formzOLD'], [r'at $v = 1$', r'at $M_{\rm current}/2$'], ['-', ':']):
            
            # ============================================================================================
            # Set up the colorbar
            # ============================================================================================
            
            TNG300_snaps = [99]
                
            root_dir = f'result/bootstrap_stats/with_{bin_type}/'
            TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
            TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])
        
            # MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
    
            # ============================================================================================
            # Plot
            # ============================================================================================
                
            # Plot TNG300       
            for snap in TNG300_snaps:
                TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                                f'med_{bin_type}', feature)
                plot_feature(simu, TNG300_z, TNG300_bin_data, TNG300_feat, [axs, cmap, norm],
                             label=label, ls=ls)
                
            # # Plot MTNG
            # for snap in MTNG_snaps:
            #     MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
            #                                               f'med_{bin_type}', feature)            
            #     plot_feature(simu, MTNG_z, MTNG_bin_data, MTNG_feat, [axs, cmap, norm])


        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xlabel(r'$z_{\rm form}$')
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        axs.legend(loc='lower right')
        
        save_dir = f'result/paper_plots/fig5/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig5_{simu}_{feature}')
        plt.close()