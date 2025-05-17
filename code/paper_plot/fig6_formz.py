import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, get_bins, plot_feature, fitting

bin_type = 'formzSub'
if bin_type == 'formz':
    xlabel = r'$z_{\rm form}$'
elif bin_type == 'formzOLD':
    xlabel = r'$z_{\rm form}$ (half mass)'
elif bin_type == 'formzSub':
    xlabel = r'$z_{\rm form}$ (half subhalo mass)'
    
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}_perMassCut/'
simus = ['Hydro']
features = ['depth','width_dimless','DWratio']


def depth(inputs, a, b):
    x, mass = inputs
    return a*x + b

def width(inputs, a, b):
    x, mass = inputs
    return a*np.log(x) + b
    
def DW(inputs, a, b, c):
    x, mass = inputs
    if bin_type == 'formzOLD':
       return (-0.3247421*x+3.42859174)/(0.62699963*np.log(x)+1.78740844) # formzOLD
    elif bin_type == 'formzSub':
        return (-0.26160318*x+3.28166956)/(0.783429*np.log(x)+1.91953979) # formzSub


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
        
        # Use mass cuts
        min_bin, max_bin, bin_width = 13, 15, 0.25
        num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

        cmap = plt.get_cmap('plasma', num_bins)
        bound = np.logspace(min_bin, max_bin+0.01, num_bins+1) 
        norm = BoundaryNorm(bound, cmap.N)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                          ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        
        # Reduce colormap ticks sf
        from matplotlib.ticker import FuncFormatter
        def custom_format(x, pos):
            return f'{x:.1f}'  
        cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
        cb.ax.tick_params(axis='x', rotation=0) 
        
        cb.set_label(r'$M_{\odot}$')
        cb._set_scale('log')

        # ============================================================================================
        # Plot
        # ============================================================================================
        
        all_x2, all_x_med, all_x_min, all_x_max = [], [], [], []
        all_y_med, all_y_min, all_y_max = [], [], []
        
        # Plot MTNG
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        MTNG_list = os.listdir(MTNG_dir)
        for fname in MTNG_list:
            MTNG_z, MTNG_x, MTNG_feat = load_stats(MTNG_dir, fname, f'med_{bin_type}', feature)    
            bin_val = 10**(10+float(fname.split('_')[3])/10)
            plot_feature(simu, bin_val, MTNG_x, MTNG_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(MTNG_x['median'])
            all_x_min.append(MTNG_x['min'])
            all_x_max.append(MTNG_x['max'])
            all_y_med.append(MTNG_feat['median'])
            all_y_min.append(MTNG_feat['min'])
            all_y_max.append(MTNG_feat['max'])
            all_x2.append(np.repeat(bin_val, len(MTNG_feat['median'])))
            
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
        
        all_x2 = np.concatenate(all_x2)
        
        # Remove NaN values
        valid_indices = ~np.isnan(all_y_med)
        all_x_med = all_x_med[valid_indices]
        all_x_min = all_x_min[valid_indices]
        all_x_max = all_x_max[valid_indices]
        all_y_med = all_y_med[valid_indices]
        all_y_min = all_y_min[valid_indices]
        all_y_max = all_y_max[valid_indices]
        all_x2 = all_x2[valid_indices]
        
        # Fit the data
        if fit_func is not None:
            popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(fit_func, 
                                                                values = [all_x_med, all_x2, all_y_med, all_y_min, all_y_max],
                                                                labels = [bin_type, 'mass', feature],
                                                                plot_info = [axs, cmap, norm], 
                                                                bootstrap=True)
            print(popt)
            print(red_chi2)
            print('')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xlabel(xlabel)
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)

        axs.legend(loc='best')
        save_dir = f'result/paper_plots/fig6_fit/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig6_{bin_type}_perMassCut_{simu}_{feature}')
        plt.close()