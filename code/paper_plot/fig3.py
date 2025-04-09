import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

bin_type = 'NFWconc'
if bin_type == 'peakHeight':
    xlabel=r'$v$'
elif bin_type == 'NFWconc':
    xlabel = r'$R_{200m}/R_s$'
elif bin_type == 'mergerz':
    xlabel = r'$z_{\rm merger}$'
    
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['width_dimless','depth']


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
            TNG300_snaps = [99]
            MTNG_snaps = [264]
            
        else:
            TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
            MTNG_snaps = [264, 237, 214, 179, 151, 129]
            
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])
       
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
        all_z = load_all_z(TNG300_dir, TNG300_snaps)
        
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('managua', len(TNG300_snaps))
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
        
        # all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
        # all_y_med, all_y_min, all_y_max = [], [], []
              
        # Plot TNG300
        if bin_type in ['mergerz', 'formz']:
            TNG300_snaps = [99]
            MTNG_snaps = [264]
        
        for snap in TNG300_snaps:
            TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                            f'med_{bin_type}', feature)
            plot_feature(simu, TNG300_z, TNG300_bin_data, TNG300_feat, [axs, cmap, norm])
            
            # # Append the data
            # all_x_med.append(TNG300_bin_data['median'])
            # all_x_min.append(TNG300_bin_data['min'])
            # all_x_max.append(TNG300_bin_data['max'])
            # all_y_med.append(TNG300_feat['median'])
            # all_y_min.append(TNG300_feat['min'])
            # all_y_max.append(TNG300_feat['max'])
            # # Duplicate z to the same length as the data
            # all_z.append([TNG300_z]*len(TNG300_bin_data['median']))
            
        # Plot MTNG
        for snap in MTNG_snaps:
            MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                      f'med_{bin_type}', feature)            
            plot_feature(simu, MTNG_z, MTNG_bin_data, MTNG_feat, [axs, cmap, norm])
            
        #     # Append the data
        #     all_x_med.append(MTNG_bin_data['median'])
        #     all_x_min.append(MTNG_bin_data['min'])
        #     all_x_max.append(MTNG_bin_data['max'])
        #     all_y_med.append(MTNG_feat['median'])
        #     all_y_min.append(MTNG_feat['min'])
        #     all_y_max.append(MTNG_feat['max'])
        #     # Duplicate z to the same length as the data
        #     all_z.append([MTNG_z]*len(MTNG_bin_data['median']))
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xlabel(xlabel)
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        # axs.legend(loc='best')
        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_{bin_type}_{simu}_{feature}')
        plt.close()