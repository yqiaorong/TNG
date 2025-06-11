""""This script plots depth (top panel) and width (bottom panel) as a function of mass at z = 0.
    Hydro simulations only."""

import os
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot the depth and width vs mass at z = 0 <<<')
print('')

root_dir = 'result/bootstrap_stats_DK14/with_mass/'
simus = ['Hydro']
features = ['width_dimless', 'depth', 'DWratio']


def mass_depth(inputs, a, b):
    mass, zval = inputs    
    return a*np.log10(mass)**b

def mass_width(inputs, a, b):
    mass, zval = inputs
    return a*np.log10(mass)**b 


for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = mass_depth
        elif feature == 'width_dimless':
            fit_func = mass_width
        
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        all_z = load_all_z(MTNG_dir, MTNG_snaps)
        all_z.append(2.1)
 
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('viridis', len(MTNG_snaps))
        bound = all_z
        norm = BoundaryNorm(bound, cmap.N)
        # cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
        #                 ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        # -----------------------------------------------------------------------------------------

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
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, 'snap_264_Rsp_stats.npy', 'med_mass', feature)
        plot_feature(simu, MTNG_z, MTNG_mass, MTNG_feat, [axs, cmap, norm])
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
        popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(fit_func, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = ['mass', 'z', feature],
                                                            plot_info = [axs, cmap, norm, '--'], 
                                                            bootstrap=True)
        print(feature)
        print(popt, red_chi2)
        print(perr)
        print('')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
            
        axs.set_xscale('log')
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        axs.set_xlabel(r'$M_{200m} /M_\odot$')
        axs.legend(loc='best')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig2_{feature}')
        plt.close()