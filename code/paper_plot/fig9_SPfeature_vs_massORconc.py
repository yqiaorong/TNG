"""This script plots the splashback features as a function of mass / concentration per 
pH cut."""

import os
import numpy as np
from func import load_stats
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
from scipy.stats import pearsonr
import pandas as pd
import seaborn as sns
plt.style.use('code/style.mplstyle')

x_type = 'mass'
if x_type == 'mass':
    x_label = r'$M_{200m} / M_{\odot}$'
else:
    x_label = r'$c$'
bin_type = f'peakHeight_per{x_type}Cut'
bin_label = r'$v$'
    
print('')
print(f'>>> Plot depth and width vs {x_type} per {bin_type} cut <<<')
print('')

root_dir = f'result/bootstrap_stats_DK14/with_{bin_type}/'
simus = ['Hydro']
features = [#'width_dimless', 
           'depth',
          # 'DWratio'
            ]


for simu in simus:
    for feature in features:

        # ============================================================================================
        # Set up the plot
        # ============================================================================================

        fig, axs = plt.subplots(1, 1, figsize=(4, 4), dpi=500, constrained_layout=False)                  

        # ============================================================================================
        # Set up the colorbar
        # ============================================================================================
        
        min_bin, max_bin, bin_width = 1, 4, 0.4 # customize
        all_bins = np.linspace(min_bin, max_bin, int((max_bin-min_bin)/bin_width+1))
        num_bins = len(all_bins) - 1

        # Set up the colorbar
        cmap = plt.get_cmap('plasma', len(all_bins))
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
        
        cb.set_label(bin_label)
        
        # ============================================================================================
        # Plot
        # ============================================================================================
                 
        # Plot MTNG
        MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'  
        MTNG_list = os.listdir(MTNG_dir)
        all_x, all_y, all_bin = [], [], []
        for fname in MTNG_list:
            print(fname)
            _, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, fname, f'med_peakHeight', feature)
            MTNG_x_val = float(fname.split('_')[1])
            
            # Duplicate z to the length of the bin data
            MTNG_x = np.repeat(MTNG_x_val, len(MTNG_bin_data['median']))
            if x_type == 'mass':
                MTNG_x = [10**(10+item) for item in MTNG_x]

            # Plot
            for i in range(len(MTNG_x)):
                axs.errorbar(MTNG_x[i], MTNG_feat['median'][i],
                             yerr=[[MTNG_feat['median'][i]-MTNG_feat['min'][i]], 
                                   [MTNG_feat['max'][i]-MTNG_feat['median'][i]]],
                            color=cmap(norm(MTNG_bin_data['median'][i])), fmt='.')
                all_x.append(MTNG_x[i])
                all_y.append(MTNG_feat['median'][i])
                all_bin.append(MTNG_bin_data['median'][i])
        all_x = np.array(all_x)
        all_y = np.array(all_y)
        all_bin = np.array(all_bin)
        
        valid_mask = (~np.isnan(all_x)) & (~np.isnan(all_y)) & (~np.isinf(all_x)) & (~np.isinf(all_y))

        all_x = all_x[valid_mask]
        all_y = all_y[valid_mask]
        all_bin = all_bin[valid_mask]
                
        # ============================================================================================
        # Correlation analysis
        # ============================================================================================
        # for bin_val in [0, 1, 2, 3]:
        #     # Select bins > 0 and < 1
        #     mask = (all_bin > bin_val) & (all_bin < bin_val+1.1)
            
        #     df = pd.DataFrame({'x': all_x[mask], 'y': all_y[mask]})
        #     sns.regplot(data=df, x='x', y='y', ci=95, line_kws={"color": cmap(norm(bin_val + 0.5))}, 
        #                scatter_kws={"s": 10})
        
        df = pd.DataFrame({'x': all_x, 'y': all_y})
        sns.regplot(data=df, x='x', y='y', ci=95, line_kws={"color": "black"}, 
                       scatter_kws={"s": 10})
        
        r_val, p_val = pearsonr(all_x, all_y)
        print(r_val, p_val)
        label_text = fr"$r = {r_val:.2f}$" + f"  (p = {p_val:.2g})"
        axs.text(0.6, 0.95, label_text,
                transform=axs.transAxes,
                fontsize=10, ha='left', va='top')
            
        # r, pval = spearmanr(all_x, all_y)

        # # Add the correlation text to the plot
        # axs.text(0.05, 0.95, f"r = {r:.2f}\n(p = {pval:.2g})",
        #         transform=axs.transAxes, fontsize=10,
        #         verticalalignment='top')
        
        # # Fit the linear regression line
        # slope, intercept, r_value, p_value, std_err = linregress(all_x, all_y)
        # ============================================================================================
        # Save the plot
        # ============================================================================================
        
        # Final edit
        axs.set_xlabel(x_label)
        if x_type == 'mass':
            axs.set_xscale('log')
        if feature == 'depth':
            Y_label = r"$\mathcal{D}$"
        elif feature == 'width_dimless':
            Y_label = r"$\mathcal{W}$"
        elif feature == 'DWratio':
            Y_label = r"$\mathcal{D}/\mathcal{W}$"
        axs.set_ylabel(Y_label)
        # axs.legend(loc='best')
        plt.tight_layout()
        save_dir = f'result/paper_plots/fig9/MTNG-Hydro/'
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(f'{save_dir}/fig9_{bin_type}_{simu}_{feature}')
        plt.close()