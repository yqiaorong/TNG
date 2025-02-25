import os
import numpy as np
from func import *
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

print('')
print(f'>>> Plot depth and width vs z (mass) <<<')
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
    # Load data 
    # ============================================================================================
    TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/snap_99/'
    TNG300_fnames = os.listdir(TNG300_dir)
    TNG300_z = sorted([float(fname.split('_')[-1][:-4]) for fname in TNG300_fnames])
    
    MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/snap_264/'
    MTNG_fnames = os.listdir(MTNG_dir)
    MTNG_z = sorted([float(fname.split('_')[-1][:-4]) for fname in MTNG_fnames])
    
    # ============================================================================================
    # Set up the colorbar
    # ============================================================================================

    min_bin, max_bin, bin_width = 11, 15, 0.5
    num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

    cmap = plt.get_cmap('vanimo', num_bins)
    bound = np.logspace(min_bin, max_bin+0.1, num_bins) 
    norm = BoundaryNorm(bound, cmap.N)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
    cb.set_label('Mass [$M_\\odot$]')
    cb._set_scale('log')

    # ============================================================================================
    # Plot
    # ============================================================================================

    for bin_val in all_bins:
        
        # Plot TNG300
        TNG300_fz, TNG300_feat = load_stats_per_bin(TNG300_dir, 'snap_99_Rsp_stats_at_{}.npy', TNG300_z, 
                                                   'formation_z', feature, 
                                                   'mass_bins', bin_val-10)
        plot_feature_vs_z(simu, 10**bin_val, TNG300_fz, TNG300_feat, [axs, cmap, norm])
            
        # Plot MTNG
        MTNG_fz, MTNG_feat = load_stats_per_bin(MTNG_dir, 'snap_264_Rsp_stats_at_{}.npy', MTNG_z,
                                               'formation_z', feature,
                                               'mass_bins', bin_val-10)
        plot_feature_vs_z(simu, 10**bin_val, MTNG_fz, MTNG_feat, [axs, cmap, norm])
        
    # Final edit
    axs.set_xlabel(r'$z_0$')
    axs.set_ylabel(r"$\mathcal{W}$")
    axs.set_xlim(0, 4)
    # ============================================================================================
    # Save the plot
    # ============================================================================================

    save_dir = f'result/paper_plots/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    plt.savefig(f'{save_dir}/fig6_formz_{simu}_{feature}')
    plt.close()