import os
import numpy as np
from func import *
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

bin_type = 'mass'
print('')
print(f'>>> Plot depth and width vs z ({bin_type}) <<<')
print('')

def mass_depth(inputs, a, d, e, f):
    zval, mass = inputs
    return e + f*zval + d*np.log10(mass)/(zval+1) + a*np.log10(mass) 

def mass_width(inputs, a, b, c, d, A, B, C, D):
    zval, mass = inputs
    return a*(zval+1)**A + b*(zval+1)**B*np.log10(mass) + c*(zval+1)**C*np.log10(mass)**2 + d*zval + D*np.exp(-zval**2)*np.log10(mass)**3


root_dir = 'result/bootstrap_stats/with_mass/'
simus = ['Hydro',]
features = ['width_dimless' ,
           'depth'
            ]
for simu in simus:
    for feature in features:
        
        if feature == 'depth':
            fit_func = mass_depth
        elif feature == 'width_dimless':
            fit_func = mass_width

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)

        # ============================================================================================
        # Load data 
        # ============================================================================================

        TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
        TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{simu}/Nboots_1024/'

        MTNG_snaps = [264, 237, 214, 179, 151, 129]
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================

        min_bin, max_bin, bin_width = 11, 15, 0.5
        num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)

        cmap = plt.get_cmap('plasma', num_bins)
        bound = np.logspace(min_bin, max_bin+0.1, num_bins) 
        norm = BoundaryNorm(bound, cmap.N)
        cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                        ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
        cb.set_label('$M_\\odot$')
        cb._set_scale('log')

        # ============================================================================================
        # Plot
        # ============================================================================================
              
        # Plot TNG300        
        for snap in TNG300_snaps:
            TNG300_z, TNG300_bin_data, TNG300_feat = load_stats(TNG300_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                            f'med_{bin_type}', feature)

            # Duplicate z to the length of the bin data
            TNG300_z = np.repeat(TNG300_z, len(TNG300_bin_data['median']))
            # Plot
            for i in range(len(TNG300_z)):
                axs.errorbar(TNG300_z[i], TNG300_feat['median'][i],
                             yerr=[[TNG300_feat['median'][i]-TNG300_feat['min'][i]], 
                                   [TNG300_feat['max'][i]-TNG300_feat['median'][i]]],
                            color=cmap(norm(TNG300_bin_data['median'][i])), fmt='.')
                 
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
        
        all_bin, all_x_med = [], []
        all_y_med, all_y_min, all_y_max = [], [], []
        
        # for bin_val in all_bins:
            
        #     # Plot TNG300
        #     TNG300_z, TNG300_feat = load_stats_per_bin(TNG300_dir, 'snap_{}_Rsp_stats.npy', TNG300_snaps, 
        #                                                'z', feature,
        #                                                'mass_bins', bin_val-10)
        #     plot_feature_vs_z(simu, 10**bin_val, TNG300_z, TNG300_feat, [axs, cmap, norm])
            
        #     # Append the data
        #     all_x_med.append(TNG300_z)
        #     all_y_med.append(TNG300_feat['median'])
        #     all_y_min.append(TNG300_feat['min'])
        #     all_y_max.append(TNG300_feat['max'])
        #     all_bin.append([10**bin_val]*len(TNG300_z))
            
        #     # Plot MTNG
        #     MTNG_z, MTNG_feat = load_stats_per_bin(MTNG_dir, 'snap_{}_Rsp_stats.npy', MTNG_snaps,
        #                                            'z', feature,
        #                                            'mass_bins', bin_val-10)
        #     plot_feature_vs_z(simu, 10**bin_val, MTNG_z, MTNG_feat, [axs, cmap, norm])
            
        #     all_x_med.append(MTNG_z)    
        #     all_y_med.append(MTNG_feat['median'])
        #     all_y_min.append(MTNG_feat['min'])
        #     all_y_max.append(MTNG_feat['max'])
        #     all_bin.append([10**bin_val]*len(MTNG_z))
  
        # Final edit
        axs.set_xlim(0, 6)
        axs.set_xlabel('z')
        if feature == 'abs_depth':
            Y_label = r"|$\mathcal{D}$|"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        else:
            Y_label = r"$\mathcal{D}$"
        axs.set_ylabel(Y_label)
        
        # # ============================================================================================
        # # Fitting
        # # ============================================================================================
        
        # # Concatenate the data to one dimension
        # all_x_med = np.concatenate(all_x_med)
    
        # all_y_med = np.concatenate(all_y_med)
        # all_y_min = np.concatenate(all_y_min)
        # all_y_max = np.concatenate(all_y_max)

        # all_bin = np.concatenate(all_bin)

        # # Fit the data
        # popt, red_chi2, y_fit, axs = fitting(all_bin, 
        #                                      all_x_med,
        #                                       all_y_med, all_y_min, all_y_max,
        #                                       fit_func, [axs, cmap, norm])
        # print(f'{simu} {feature} popt: {popt}')
        # print(f'{simu} {feature} red_chi2: {red_chi2}')
        # axs.legend()

        # ============================================================================================
        # Save the plot
        # ============================================================================================

        save_dir = f'result/paper_plots/fig4/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig4_mass_{simu}_{feature}')
        plt.close()