import os
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot depth and width vs accretion rate <<<')
print('')

root_dir = 'result/bootstrap_stats/with_accret/'
simus = ['Hydro','DM']
features = ['width_dimless', 
            #'abs_depth',
           # 'depth'
            ]
upper_limit = 6
    

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================
def accret_depth(inputs, a, b, c, e, f):
    Gamma, zval = inputs
    
    return a*Gamma + b*Gamma**2 + c*Gamma**3 + e/(zval+1) + f/Gamma

def accret_width(inputs, a, b, c, d, e, f,):
    Gamma, zval = inputs
    # return a*zval**2*Gamma + b*zval**2 + c*zval + d*Gamma + e*Gamma*zval + f
    return  a*Gamma + b*Gamma**2 + c*Gamma**3 + d/(zval+1) + f*zval + e

for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = accret_depth
        elif feature == 'width_dimless':
            fit_func = accret_width

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'
        TNG300_z_i, TNG300_z_f = load_all_z(TNG300_dir, [min(TNG300_snaps), max(TNG300_snaps)])

        all_z = load_all_z(TNG300_dir, TNG300_snaps)

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

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

            TNG300_z, TNG300_accret, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy', 
                                                            'med_accret', feature)
            # Filter out data greater than 6
            mask = TNG300_accret['median'] <= upper_limit
            TNG300_accret['median'] = TNG300_accret['median'][mask]
            TNG300_accret['min'] = TNG300_accret['min'][mask]
            TNG300_accret['max'] = TNG300_accret['max'][mask]
            TNG300_feat['median'] = TNG300_feat['median'][mask]
            TNG300_feat['min'] = TNG300_feat['min'][mask]
            TNG300_feat['max'] = TNG300_feat['max'][mask]
            
            plot_feature(simu, TNG300_z, TNG300_accret, TNG300_feat, [axs, cmap, norm])
            
            # Append the data
            all_x_med.append(TNG300_accret['median'])
            all_x_min.append(TNG300_accret['min'])
            all_x_max.append(TNG300_accret['max'])
            all_y_med.append(TNG300_feat['median'])
            all_y_min.append(TNG300_feat['min'])
            all_y_max.append(TNG300_feat['max'])
            # Duplicate z to the same length as the data
            all_z.append([TNG300_z]*len(TNG300_accret['median']))

        # # Plot MTNG
        # for snap in MTNG_snaps:

        #     MTNG_z, MTNG_accret, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
        #                                               'med_accret', feature)
        #     # Filter out data greater than 6
        #     mask = MTNG_accret['median'] <= upper_limit
        #     MTNG_accret['median'] = MTNG_accret['median'][mask]
        #     MTNG_accret['min'] = MTNG_accret['min'][mask]
        #     MTNG_accret['max'] = MTNG_accret['max'][mask]
        #     MTNG_feat['median'] = MTNG_feat['median'][mask]
        #     MTNG_feat['min'] = MTNG_feat['min'][mask]
        #     MTNG_feat['max'] = MTNG_feat['max'][mask]
            
        #     plot_feature(simu, MTNG_z, MTNG_accret, MTNG_feat, [axs, cmap, norm])
            
        #     # Append the data
        #     all_x_med.append(MTNG_accret['median'])
        #     all_x_min.append(MTNG_accret['min'])
        #     all_x_max.append(MTNG_accret['max'])
        #     all_y_med.append(MTNG_feat['median'])
        #     all_y_min.append(MTNG_feat['min'])
        #     all_y_max.append(MTNG_feat['max'])
        #     # Duplicate z to the same length as the data
        #     all_z.append([MTNG_z]*len(MTNG_accret['median']))
        
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
        
        # Fit the data
        popt, red_chi2, y_fit, axs = fitting(all_z, 
                                              all_x_med, all_x_min, all_x_max, 
                                              all_y_med, all_y_min, all_y_max,
                                              fit_func, [axs, cmap, norm])
        print(f'{simu} {feature} popt: {popt}')
        print(f'{simu} {feature} red_chi2: {red_chi2}')
        
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        # axs.set_xlim(0, 6)
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
            axs.set_ylim(1, 4)
        axs.set_ylabel(Y_label)
        axs.set_xlabel(r'$\Gamma$')
        axs.legend()

        save_dir = f'result/paper_plots/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig3_accret_{simu}_{feature}')
        plt.close()