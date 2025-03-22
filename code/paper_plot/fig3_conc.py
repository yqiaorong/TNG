import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot depth and width vs mass <<<')
print('')

root_dir = 'result/bootstrap_stats/with_conc/'
simus = ['Hydro','DM']
features = [# abs_depth', 
            'width_dimless',
           'depth', 
            ]

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================
def conc_depth(inputs, a, d, e, f):
    x, zval = inputs    
    return e + f*zval + d*np.log10(x)/(zval+1) + a*np.log10(x) 

def conc_width(inputs, a, b, c, d, A, B, C, D):
    x, zval = inputs
    return a*(zval+1)**A + b*(zval+1)**B*np.log10(x) + c*(zval+1)**C*np.log10(x)**2 + d*zval + D*np.exp(-zval**2)*np.log10(x)**3

for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = conc_depth
        elif feature == 'width_dimless':
            fit_func = conc_width
            
        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])
        
        if simu == 'DM':
            MTNG_snaps = [264, 237, 214, 179]
        else:
            MTNG_snaps = [264, 237, 214, 179, 151]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'
        
        all_z = load_all_z(TNG300_dir, TNG300_snaps)
        
        # Set up the colorbar
        # -----------------------------------------------------------------------------------------
        cmap = plt.get_cmap('managua', len(TNG300_snaps))
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
              
        # Plot TNG300
        for snap in TNG300_snaps:
            TNG300_z, TNG300_conc, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                            'med_conc', feature)
            plot_feature(simu, TNG300_z, TNG300_conc, TNG300_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(TNG300_conc['median'])
            all_x_min.append(TNG300_conc['min'])
            all_x_max.append(TNG300_conc['max'])
            all_y_med.append(TNG300_feat['median'])
            all_y_min.append(TNG300_feat['min'])
            all_y_max.append(TNG300_feat['max'])
            # Duplicate z to the same length as the data
            all_z.append([TNG300_z]*len(TNG300_conc['median']))
            
        # Plot MTNG
        for snap in MTNG_snaps:
            MTNG_z, MTNG_conc, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                      'med_conc', feature)
            # Remove extremely small points
            remove_idx = np.where(MTNG_conc['median'] < 0.01)
            MTNG_conc = {k: np.delete(v, remove_idx) for k, v in MTNG_conc.items()}
            MTNG_feat = {k: np.delete(v, remove_idx) for k, v in MTNG_feat.items()}
            
            plot_feature(simu, np.round(MTNG_z, 1), MTNG_conc, MTNG_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(MTNG_conc['median'])
            all_x_min.append(MTNG_conc['min'])
            all_x_max.append(MTNG_conc['max'])
            all_y_med.append(MTNG_feat['median'])
            all_y_min.append(MTNG_feat['min'])
            all_y_max.append(MTNG_feat['max'])
            # Duplicate z to the same length as the data
            all_z.append([MTNG_z]*len(MTNG_conc['median']))

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
        
        # # Fit the data
        # popt, red_chi2, y_fit, axs = fitting(all_z, 
        #                                       all_x_med, 
        #                                       all_y_med, all_y_min, all_y_max,
        #                                       fit_func, [axs, cmap, norm],
        #                                       all_x_min, all_x_max, )
        # print(popt)
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        axs.set_xscale('log')
        axs.set_xlabel(r"$R_s/R_{200m}$")
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        axs.legend(loc='best')
        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_conc_{simu}_{feature}')
        plt.close()