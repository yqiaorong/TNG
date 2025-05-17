import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, plot_feature, fitting

bin_type = 'formzSub'
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['width_dimless', 'depth', 'DWratio']

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================

def depth(inputs, a, b):
    x, zval = inputs
    return a*x + b

def width(inputs, a, b):
    x, zval = inputs
    return a*np.log(x) + b
    
def DW(inputs, a, b, c):
    x, zval = inputs
    if bin_type == 'formzOLD':
        return (-0.34200887*x+3.14556083)/(0.5975754*np.log(x)+2.04735604) # formzOLD
    elif bin_type == 'formzSub':
        return (-0.46966875*x+3.19050738)/(0.5514771*np.log(x)+2.03946774) # formzSub


for simu in simus:
    for feature in features:
        print(feature)
        if feature == 'depth':
            fit_func = depth
        elif feature == 'width_dimless':
            fit_func = width
        elif feature == 'DWratio':
            fit_func = DW
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.3), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        MTNG_snaps = [264, 
                    #  237, 214, 179, 151, 129
                      ]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
        # all_z = load_all_z(MTNG_dir, MTNG_snaps)
        # all_z.append(2.1)
        all_z = np.array([0, 1, 2])
        
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('viridis', len(all_z))
        norm = BoundaryNorm(all_z, cmap.N)
        # cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
        #                 ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        # # -----------------------------------------------------------------------------------------

        # # Reduce colormap ticks sf
        # from matplotlib.ticker import FuncFormatter
        # def custom_format(x, pos):
        #     return f'{x:.1f}'  
        # cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        # cb.ax.tick_params(axis='x', rotation=70) 
        # cb.set_label('z')

        # ============================================================================================
        # Plot
        # ============================================================================================
        
        all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
        all_y_med, all_y_min, all_y_max = [], [], []
            
        # Plot MTNG
        for snap in MTNG_snaps:
            MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                      f'med_{bin_type}', feature)
            plot_feature(simu, np.round(MTNG_z, 1), MTNG_bin_data, MTNG_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(MTNG_bin_data['median'])
            all_x_min.append(MTNG_bin_data['min'])
            all_x_max.append(MTNG_bin_data['max'])
            all_y_med.append(MTNG_feat['median'])
            all_y_min.append(MTNG_feat['min'])
            all_y_max.append(MTNG_feat['max'])
            # Duplicate z to the same length as the data
            all_z.append([MTNG_z]*len(MTNG_bin_data['median']))
        
        # ============================================================================================
        # Fitting
        # ============================================================================================
        
        # Concatenate the data to one dimension
        all_x_med = np.concatenate(all_x_med)

        all_x_min = np.concatenate(all_x_min)
        all_x_max = np.concatenate(all_x_max)
            
        all_y_med = np.concatenate(all_y_med)
        all_y_min = np.concatenate(all_y_min)
        all_y_max = np.concatenate(all_y_max)
        
        all_z = np.concatenate(all_z)
        
        # Remove NaN values
        valid_indices = ~np.isnan(all_y_med)
        all_x_med = all_x_med[valid_indices]
        all_x_min = all_x_min[valid_indices]
        all_x_max = all_x_max[valid_indices]
        all_y_med = all_y_med[valid_indices]
        all_y_min = all_y_min[valid_indices]
        all_y_max = all_y_max[valid_indices]
        all_z = all_z[valid_indices]
        
        # Fit the data
        if fit_func is not None:
            popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(fit_func, 
                                                                values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                                labels = [bin_type, 'z', feature],
                                                                plot_info = [axs, cmap, norm], 
                                                                bootstrap=True)
            print(popt)
            print(perr)
            print(red_chi2)
            
            # Save fitted data
            fitted_data_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
            if not os.path.exists(fitted_data_dir):
                os.makedirs(fitted_data_dir)
            np.save(fitted_data_dir + f'{bin_type}_{feature}_fitted_data.npy', fitted_data)
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xlabel(r"$\nu$")
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
            if bin_type == 'NFWconc':
                axs.set_ylim(1, 3.8)
        axs.set_ylabel(Y_label)
        axs.legend(loc='best')
        
        save_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_{bin_type}_{simu}_{feature}')
        plt.close()