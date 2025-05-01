import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, load_all_z, plot_feature, fitting

bin_type = 'NFWconc'
print('')
print(f'>>> Plot depth and width vs {bin_type} <<<')
print('')

root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simus = ['Hydro']
features = ['width_dimless', 'depth', 'DWratio']

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================


def conc_depth(inputs, a, b, c):
    conc, zval = inputs    
    return a*conc**b / (zval + 1)**c

def conc_width(inputs, a,b, d):
    conc, zval = inputs
    return a*conc*np.exp(b*(conc-d))#*(zval + 1)**c
   # return a* np.exp(c*np.log10(mass)) / (zval + 1)**d
    
def conc_DW(inputs, a, b, c):
    conc, zval = inputs
    return a*conc**b / (zval + 1)**c


for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = conc_depth
        elif feature == 'width_dimless':
            fit_func = conc_width
        elif feature == 'DWratio':
            fit_func = conc_DW
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        # TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
        # TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        # TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
        all_z = load_all_z(MTNG_dir, MTNG_snaps)
        
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('viridis', len(MTNG_snaps))
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
        
        all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
        all_y_med, all_y_min, all_y_max = [], [], []
              
        # # Plot TNG300
        # for snap in TNG300_snaps:
        #     TNG300_z, TNG300_mass, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
        #                                                    f'med_{bin_type}', feature)
        #     plot_feature(simu, TNG300_z, TNG300_mass, TNG300_feat, [axs, cmap, norm])
            
        #     # Append the data
        #     all_x_med.append(TNG300_mass['median'])
        #     all_x_min.append(TNG300_mass['min'])
        #     all_x_max.append(TNG300_mass['max'])
        #     all_y_med.append(TNG300_feat['median'])
        #     all_y_min.append(TNG300_feat['min'])
        #     all_y_max.append(TNG300_feat['max'])
        #     # Duplicate z to the same length as the data
        #     all_z.append([TNG300_z]*len(TNG300_mass['median']))
            
        # Plot MTNG
        for snap in MTNG_snaps:
            MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                      f'med_{bin_type}', feature)
            plot_feature(simu, np.round(MTNG_z, 1), MTNG_mass, MTNG_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(MTNG_mass['median'])
            all_x_min.append(MTNG_mass['min'])
            all_x_max.append(MTNG_mass['max'])
            all_y_med.append(MTNG_feat['median'])
            all_y_min.append(MTNG_feat['min'])
            all_y_max.append(MTNG_feat['max'])
            # Duplicate z to the same length as the data
            all_z.append([MTNG_z]*len(MTNG_mass['median']))
        
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
        popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(bin_type, feature, 
                                                                all_z, 
                                                                all_x_med, 
                                                                all_y_med, all_y_min, all_y_max,
                                                                fit_func, [axs, cmap, norm],
                                                                all_x_min, all_x_max, )
        print(popt, red_chi2)
        
        
        # Save fitted data
        fitted_data_dir = f'result/paper_plots/fig4/MTNG-Hydro/'
        if not os.path.exists(fitted_data_dir):
            os.makedirs(fitted_data_dir)
        np.save(fitted_data_dir + f'{bin_type}_{feature}_fitted_data.npy', fitted_data)
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
    
        axs.set_xlabel("c")
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        axs.legend(loc='best')
        save_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_{bin_type}_{simu}_{feature}')
        plt.close()