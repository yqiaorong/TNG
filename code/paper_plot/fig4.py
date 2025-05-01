import os
import numpy as np
from func import get_bins, load_stats
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

bin_type = 'accretions'
if bin_type == 'peakHeight':
    label=r'$v$'
elif bin_type == 'NFWconc':
    label = r'$c$'
elif bin_type == 'mergerz':
    label = r'$z_{\rm merger}$'
elif bin_type == 'formz':
    label = r'$z_{\rm form}$'
elif bin_type == 'formzOLD':
    label = r'$z_{\rm form}$ (half mass)'
elif bin_type in ['accretions', 'accretionsOLD', 'accret']:
    label = r'$\Gamma$'
elif bin_type == 'mass':
    label = r'$M_{\odot}$'
    

print('')
print(f'>>> Plot depth and width vs z ({bin_type}) <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['width_dimless' , 'depth', 'DWratio']
for simu in simus:
    for feature in features:

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Load data 
        # ============================================================================================
        
        if bin_type == 'formz':
            TNG300_snaps = [99, 78, 67, 50]
        elif bin_type == 'formzOLD':
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25]
        elif bin_type == 'accretions':
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13]
        else:
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'            

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================
        
        if bin_type == 'peakHeight':
            min_bin, max_bin, bin_width = 0, 4, 0.5 # customize
        elif bin_type == 'mass':
            min_bin, max_bin, bin_width = 13, 15.5, 0.5
        elif bin_type == 'accretions':
            min_bin, max_bin, bin_width = 0, 6, 0.5
        elif bin_type == 'NFWconc':
            min_bin, max_bin, bin_width = 0, 18, 2
        elif bin_type == 'formz':
            min_bin, max_bin, bin_width = 0, 2, 0.2
        elif bin_type == 'formzOLD':
            min_bin, max_bin, bin_width = 0, 3.5, 0.4
        num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

        cmap = plt.get_cmap('plasma', num_bins)
        if bin_type == 'mass':
            bound = np.logspace(min_bin, max_bin+0.01, num_bins+1) 
        else:
            bound = np.linspace(min_bin, max_bin+0.01, num_bins+1) 
        norm = BoundaryNorm(bound, cmap.N)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                          ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        
        # Reduce colormap ticks sf
        from matplotlib.ticker import FuncFormatter
        def custom_format(x, pos):
            return f'{x:.1f}'  
        cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        cb.ax.tick_params(axis='x', rotation=0) 
        
        cb.set_label(label)
        if bin_type == 'mass':
            cb._set_scale('log')

        # ============================================================================================
        # Plot
        # ============================================================================================
              
        # # Plot TNG300        
        # for snap in TNG300_snaps:
        #     TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
        #                                                     f'med_{bin_type}', feature)

        #     # Duplicate z to the length of the bin data
        #     TNG300_z = np.repeat(TNG300_z, len(TNG300_bin_data['median']))
        #     # Plot
        #     for i in range(len(TNG300_z)):
        #         axs.errorbar(TNG300_z[i], TNG300_feat['median'][i],
        #                      yerr=[[TNG300_feat['median'][i]-TNG300_feat['min'][i]], 
        #                            [TNG300_feat['max'][i]-TNG300_feat['median'][i]]],
        #                     color=cmap(norm(TNG300_bin_data['median'][i])), fmt='.')
                 
        # Plot MTNG
        for snap in MTNG_snaps:
            MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                            f'med_{bin_type}', feature)

            # Duplicate z to the length of the bin data
            MTNG_z = np.repeat(MTNG_z, len(MTNG_bin_data['median']))
            # Plot
            for i in range(len(MTNG_z)):
                axs.errorbar(MTNG_z[i], MTNG_feat['median'][i],
                             yerr=[[MTNG_feat['median'][i]-MTNG_feat['min'][i]], 
                                   [MTNG_feat['max'][i]-MTNG_feat['median'][i]]],
                            color=cmap(norm(MTNG_bin_data['median'][i])), fmt='.')
        
        # ============================================================================================
        # Load fitted data and plot
        # ============================================================================================
        # fitted_dir = f'result/paper_plots/fig4/MTNG-Hydro/'
        # fitted_data = np.load(fitted_dir + f'{bin_type}_{feature}_fitted_data.npy', allow_pickle=True).item()
        
        # for bin_start, bin_end in zip(all_bins[:-1], all_bins[1:]):
        #     if bin_type == 'mass':
        #         select_idx = np.where((np.log10(fitted_data[bin_type]) >= bin_start) &
        #                                 (np.log10(fitted_data[bin_type]) < bin_end))[0]
        #         color_val = np.mean(np.log10(fitted_data[bin_type][select_idx]))
        #         plot_x = fitted_data['z'][select_idx]
        #         plot_y = fitted_data[feature][select_idx]
        #         axs.plot(plot_x, plot_y, color=cmap(norm(10**color_val)), ls='--') 
        #     else:
        #         select_idx = np.where((fitted_data[bin_type] >= bin_start) &
        #                               (fitted_data[bin_type] < bin_end))[0]
        #         if select_idx.size == 0:
        #             pass
        #         else:
        #             color_val = np.mean(fitted_data[bin_type][select_idx])
        #             plot_x = fitted_data['z'][select_idx]
        #             plot_y = fitted_data[feature][select_idx]
        #             axs.plot(plot_x, plot_y, color=cmap(norm(color_val)), ls='--') 
        
        
        # Final edit
        axs.set_xlabel('z')
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        axs.legend(loc='best')
        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/fig4/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig4_{bin_type}_{simu}_{feature}')
        plt.close()