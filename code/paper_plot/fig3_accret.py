import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot depth and width vs accretion rate <<<')
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
        # Set up the colorbar
        # ============================================================================================

        TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])

        all_z = load_all_z(TNG300_dir, TNG300_snaps)

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('managua', len(TNG300_snaps))
        bound = all_z
        norm = BoundaryNorm(bound, cmap.N)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                        ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        # -----------------------------------------------------------------------------------------

        # Reduce colormap ticks sf
        from matplotlib.ticker import FuncFormatter
        def custom_format(x, pos):
            return f'{x:.1f}'  
        cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        cb.ax.tick_params(axis='x', rotation=70) 
        cb.set_label('z')

        # ============================================================================================
        # Plot
        # ============================================================================================

        # Plot TNG300
        for snap in TNG300_snaps:  

            TNG300_z, TNG300_mass, TNG300_feat = load_stats(TNG300_dir, snap, 'med_accret', feature)
            plot_feature(simu, TNG300_z, TNG300_mass, TNG300_feat, [axs, cmap, norm])

        # Plot MTNG
        for snap in MTNG_snaps:

            MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, snap, 'med_accret', feature)
            plot_feature(simu, MTNG_z, MTNG_mass, MTNG_feat, [axs, cmap, norm])

        axs.set_xlim(0, 6)
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        axs.set_xlabel(r'$\Gamma$')

        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_accret_{simu}_{feature}')
        plt.close()