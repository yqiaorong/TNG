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

root_dir = 'result/bootstrap_stats/with_accret/'
simus = ['Hydro', 'DM']
features = ['width_dimless' , 'abs_depth', 'depth']
for simu in simus:
    for feature in features:

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

        min_bin, max_bin, bin_width = 0, 6, 1
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

        for bin_val in all_bins:
            
            # Plot TNG300
            TNG300_z, TNG300_feat = load_stats_per_bin(TNG300_dir, 'snap_{}_Rsp_stats.npy', TNG300_snaps, 
                                                       'z', feature,
                                                       'accret_bins', bin_val)
            plot_feature_vs_z(simu, bin_val, TNG300_z, TNG300_feat, [axs, cmap, norm])
                
            # Plot MTNG
            MTNG_z, MTNG_feat = load_stats_per_bin(MTNG_dir, 'snap_{}_Rsp_stats.npy', MTNG_snaps,
                                                   'z', feature,
                                                   'accret_bins', bin_val)
            plot_feature_vs_z(simu, bin_val, MTNG_z, MTNG_feat, [axs, cmap, norm])
            
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
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig4_accret_{simu}_{feature}')
        plt.close()