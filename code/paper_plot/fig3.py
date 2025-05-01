import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, load_all_z, plot_feature

bin_type = 'NFWconc'
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
elif bin_type in ['accretions', 'accretionsOLD', 'accret']:
    xlabel = r'$\Gamma$'
elif bin_type == 'mass':
    xlabel = r'$M_{\odot}$'
    
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['width_dimless','depth', 'DWratio']


for simu in simus:
    for feature in features:
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================
        
        if bin_type in ['mergerz', 'formz']:
            TNG300_snaps = [99, 78, 67, 50]
            MTNG_snaps = [264]
        elif bin_type in ['formzOLD', 'formzSub']:
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25]
            MTNG_snaps = [264]
        elif bin_type in ['accretions']:
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
            MTNG_snaps = [264, 237, 214, 179, 151, 129]
        else:
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
            MTNG_snaps = [264, 
                        237, 214, 179, 
                          151, 
                         129
                          ]
            
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
       
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
        all_z = load_all_z(MTNG_dir, MTNG_snaps)
        all_z.append(2.1)
        # if bin_type in ['mergerz', 'formz']:
        #     all_z = np.linspace(min(all_z), max(all_z), 10)
            
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('viridis', len(all_z))
        if len(all_z) == 1:
            bound = [all_z[0], all_z[0]+0.1]
        else:
            bound = all_z
        norm = BoundaryNorm(bound, cmap.N)
        
        if len(all_z) > 1:
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
              
        # # Plot TNG300       
        # for snap in TNG300_snaps:
        #     TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
        #                                                     f'med_{bin_type}', feature)
        #     plot_feature(simu, TNG300_z, TNG300_bin_data, TNG300_feat, [axs, cmap, norm])
            
        # Plot MTNG
        for snap in MTNG_snaps:
            print(snap)
            MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                      f'med_{bin_type}', feature)            
            plot_feature(simu, MTNG_z, MTNG_bin_data, MTNG_feat, [axs, cmap, norm], label=f'z={MTNG_z:.1f}')

        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        if feature == 'DWratio':
            axs.set_ylim(0, 6)
        if bin_type == 'mass':
            axs.set_xscale('log')
        axs.set_xlabel(xlabel)
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        # axs.legend(loc='best')
        
        save_dir = f'result/paper_plots/fig3/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_{bin_type}_{simu}_{feature}')
        plt.close()