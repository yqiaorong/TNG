import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot depth and width vs accretion rate <<<')
print('')

root_dir = 'result/bootstrap_stats/with_formation_z/'
simus = ['Hydro', 'DM']
feature = 'width_dimless' 

for simu in simus:
    # ============================================================================================
    # Set up the plot
    # ============================================================================================

    fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

    # ============================================================================================
    # Set up the colorbar
    # ============================================================================================

    TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/snap_99/'
    TNG300_fnames = os.listdir(TNG300_dir)
    TNG300_z = sorted([float(fname.split('_')[-1][:-4]) for fname in TNG300_fnames])

    MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/snap_264/'
    MTNG_fnames = os.listdir(MTNG_dir)
    MTNG_z = sorted([float(fname.split('_')[-1][:-4]) for fname in MTNG_fnames])
    
    # Set up the colorbar
    # -----------------------------------------------------------------------------------------
    cmap = plt.get_cmap('managua', len(TNG300_z))
    bound = TNG300_z
    norm = BoundaryNorm(bound, cmap.N)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                      ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
    # Reduce colormap ticks sf
    from matplotlib.ticker import FuncFormatter
    def custom_format(x, pos):
        return f'{x:.1f}'  
    cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
    cb.ax.tick_params(axis='x', rotation=70) 
    cb.set_label('z')
    # -----------------------------------------------------------------------------------------
    
    # ============================================================================================
    # Plot
    # ============================================================================================

    # Plot TNG300
    for z in TNG300_z:  
        TNG300_z, TNG300_mass, TNG300_feat = load_stats(TNG300_dir, f'snap_99_Rsp_stats_at_{z}.npy', 
                                                        'med_mass', feature, 'formation_z')
        plot_feature(simu, TNG300_z, TNG300_mass, TNG300_feat, [axs, cmap, norm])

    # Plot MTNG
    for z in MTNG_z:
        MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, f'snap_264_Rsp_stats_at_{z}.npy', 
                                                  'med_mass', feature, 'formation_z')
        plot_feature(simu, MTNG_z, MTNG_mass, MTNG_feat, [axs, cmap, norm])

    axs.set_ylabel(r"$\mathcal{W}$")
    axs.set_xlabel(r"$M_\odot$")
    axs.set_xscale('log')
    
    # ============================================================================================
    # Save the plot
    # ============================================================================================

    save_dir = f'result/paper_plots/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.savefig(f'{save_dir}/fig5_formz_{simu}_{feature}')
    plt.close()