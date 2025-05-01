""""This script plots depth (top panel) and width (bottom panel) as a function of mass at z = 0.
    Hydro simulations only."""

import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot the absolute depth and width vs mass at z = 0 <<<')
print('')

root_dir = 'result/bootstrap_stats/with_mass/'
simus = ['Hydro']
features = ['width_dimless', 'depth', 'DWratio']




for simu in simus:
    for feature in features:
        
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        # TNG300_snaps = [99, 78, 67, 50, 40, 33, 21, 17, 8]
        # TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simus[0]}/Nboots_1024/'
        # TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])


        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        all_z = load_all_z(MTNG_dir, MTNG_snaps)

        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('viridis', len(MTNG_snaps))
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
        
        # # Plot TNG300
        # TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        # TNG300_z, TNG300_mass, TNG300_depths = load_stats(TNG300_dir, 'snap_99_Rsp_stats.npy', 'med_mass', 'depth')
        # _, _, TNG300_widths = load_stats(TNG300_dir, 'snap_99_Rsp_stats.npy', 'med_mass', 'width_dimless')
        # plot_feature(simu, TNG300_z, TNG300_mass, TNG300_depths, [axs[0], cmap, norm])
        # plot_feature(simu, TNG300_z, TNG300_mass, TNG300_widths, [axs[1], cmap, norm])
                
        # Plot MTNG
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        MTNG_z, MTNG_mass, MTNG_feature = load_stats(MTNG_dir, 'snap_264_Rsp_stats.npy', 'med_mass', feature)
        plot_feature(simu, MTNG_z, MTNG_mass, MTNG_feature, [axs, cmap, norm])

        axs.set_xscale('log')
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        axs.set_xlabel('$M_\\odot$')

        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig2_{feature}')
        plt.close()