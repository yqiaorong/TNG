import os
import numpy as np
from func import get_bins, load_stats, fitting
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

bin_type = 'accretions'
label = r'$\Gamma$'
    
print('')
print(f'>>> Plot depth and width vs z ({bin_type}) <<<')
print('')


def accret_depth(inputs, a, b, c):
    zval, accret = inputs    
    return a*zval**b + c

def accret_width(inputs, a, b, c):
    zval, accret = inputs   
    return a*zval**b + c
    
def accret_DW(inputs, a, b, c):
    zval, accret = inputs   
    return a*zval**b + c


root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['depth', 'width_dimless', 'DWratio']
for simu in simus:
    for feature in features:
        print(feature)
        
        if feature == 'depth':
            fit_func = accret_depth
        elif feature == 'width_dimless':
            fit_func = accret_width
        elif feature == 'DWratio':
            fit_func = accret_DW

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Load data 
        # ============================================================================================

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'            

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================
        
        min_bin, max_bin, bin_width = 0, 6, 0.5
        num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

        cmap = plt.get_cmap('plasma', num_bins)
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

        # ============================================================================================
        # Plot
        # ============================================================================================
        
        all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
        all_y_med, all_y_min, all_y_max = [], [], []
        
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
                # Append the data
                all_z.append(MTNG_z[i])
                all_x_med.append(MTNG_bin_data['median'][i])
                all_x_min.append(MTNG_bin_data['min'][i])
                all_x_max.append(MTNG_bin_data['max'][i])
                all_y_med.append(MTNG_feat['median'][i])
                all_y_min.append(MTNG_feat['min'][i])
                all_y_max.append(MTNG_feat['max'][i])
                
        # ============================================================================================
        # Fitting
        # ============================================================================================
        
        # Concatenate the data to one dimension
        all_x_med = np.asarray(all_x_med)

        all_x_min = np.asarray(all_x_min)
        all_x_max = np.asarray(all_x_max)
            
        all_y_med = np.asarray(all_y_med)
        all_y_min = np.asarray(all_y_min)
        all_y_max = np.asarray(all_y_max)
        
        all_z = np.asarray(all_z)
        
        # Remove nan
        mask = ~np.isnan(all_x_med)
        all_x_med = all_x_med[mask]
        all_x_min = all_x_min[mask]
        all_x_max = all_x_max[mask]
        all_y_med = all_y_med[mask]
        all_y_min = all_y_min[mask]
        all_y_max = all_y_max[mask]
        all_z = all_z[mask]
        
        # Fit the data
        popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(fit_func, 
                                                                values = [all_z, all_x_med, all_y_med, all_y_min, all_y_max],
                                                                labels = ['z', bin_type, feature],
                                                                plot_info = [axs, cmap, norm], 
                                                                bootstrap=True)
        print(popt, red_chi2)
        print('')
        
        # # ============================================================================================
        # # Load fitted data and plot
        # # ============================================================================================
        
        # fitted_dir = f'result/paper_plots/fig4_fit/MTNG-Hydro/'
        # fitted_data = np.load(fitted_dir + f'{bin_type}_{feature}_fitted_data.npy', allow_pickle=True).item()
        # fitted_bins = np.unique(fitted_data[bin_type])
        # print(fitted_data)
        # print(fitted_bins)
        # for bin_start, bin_end in zip(fitted_bins[:-1], fitted_bins[1:]):
        #     select_idx = np.where((fitted_data[bin_type] >= bin_start) & (fitted_data[bin_type] < bin_end))[0]
        #     if select_idx.size == 0:
        #         pass
        #     else:
        #         color_val = np.mean(fitted_data[bin_type][select_idx])
        #         plot_x = fitted_data['z'][select_idx]
        #         plot_y = fitted_data[feature][select_idx]
        #         axs.plot(plot_x, plot_y, color=cmap(norm(color_val)), ls='--') 
        
        
        # Final edit
        axs.set_xlabel('z')
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
            # if bin_type == 'peakHeight':
            #     axs.set_ylim(2.3, 4)
            # elif bin_type == 'NFWconc':
            #     axs.set_ylim(2.3, 3.5)
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
            # if bin_type == 'NFWconc':
            #     axs.set_ylim(1, 3.8)
        axs.set_ylabel(Y_label)
        axs.legend(loc='best')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/fig4_fit/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig4_{bin_type}_{simu}_{feature}')
        plt.close()