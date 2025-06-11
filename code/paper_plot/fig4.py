import os
import numpy as np
from func import get_bins, load_stats
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FormatStrFormatter
plt.style.use('code/style.mplstyle')

bin_type = 'mass'
print('')
print(f'>>> Plot depth and width vs z ({bin_type}) <<<')
print('')



if bin_type == 'peakHeight':
    label=r'$v$'
    min_bin, max_bin, bin_width = 0, 4, 0.5 # customize
elif bin_type == 'NFWconc':
    label = r'$c$'
    min_bin, max_bin, bin_width = 0, 18, 2
# elif bin_type == 'mergerz':
#     label = r'$z_{\rm merger}$'
elif bin_type == 'formz':
    label = r'$z_{\rm form}$'
    min_bin, max_bin, bin_width = 0, 2, 0.2
elif bin_type == 'formzOLD':
    label = r'$z_{\rm form}$ (half mass)'
    min_bin, max_bin, bin_width = 0, 3.5, 0.4
# elif bin_type in ['accretions', 'accretionsOLD', 'accret']:
#     label = r'$\Gamma$'
#     min_bin, max_bin, bin_width = 0, 6, 0.5
elif bin_type == 'mass':
    label = r"$M_{200m} / M_\odot$"
    min_bin, max_bin, bin_width = 13, 15.5, 0.5
num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)
    


root_dir = f'result/bootstrap_stats_DK14/with_{bin_type}/'
simu = 'Hydro'
features = ['depth', 'width_dimless', #'DWratio'
            ]
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$", # r"$\mathcal{D}/\mathcal{W}$"
           ]

fig, axs = plt.subplots(2, 1, figsize=(4, 6), dpi=500, sharex=True, constrained_layout=True) 


# Set up the colorbar
cmap = plt.get_cmap('plasma', len(all_bins))
if bin_type == 'mass':
    bound = np.logspace(min(all_bins), max(all_bins)+0.01, len(all_bins)+1) 
else:
    bound = np.linspace(min(all_bins), max(all_bins)+0.01, len(all_bins)+1) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=axs[1], orientation='horizontal', spacing='proportional', ticks=bound)

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.ax.tick_params(axis='x', rotation=0) 

cb.set_label(label)
if bin_type == 'mass':
    cb._set_scale('log')
        
        
# ============================================================================================
# Load data 
# ============================================================================================

MTNG_snaps = [264, 237, 214, 179, 151, 129]
MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'      

for ifeat, feature in enumerate(features):
    
    if bin_type == 'NFWconc' and feature == 'depth':
        pass
    else:
        # ============================================================================================
        # Load fitted data and plot
        # ============================================================================================
        
        fitted_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
        fitted_data = np.load(fitted_dir + f'{bin_type}_{feature}_fitted_data.npy', allow_pickle=True).item()
        
        fitted_vals_bins = np.unique(fitted_data[bin_type])
        fitted_vals      = np.concatenate(fitted_data[bin_type])
        fitted_z         = np.concatenate(fitted_data['z'])
        fitted_y         = np.concatenate(fitted_data[feature])
        del fitted_data

        # Plot the fitted data
        for bin_start, bin_end in zip(fitted_vals_bins[:-1], fitted_vals_bins[1:]):
            select_idx = np.where((fitted_vals >= bin_start) & (fitted_vals < bin_end))[0]
            if select_idx.size == 0:
                pass
            else:
                color_val = np.mean(fitted_vals[select_idx])

                plot_x = fitted_z[select_idx]
                plot_y = fitted_y[select_idx]
                axs[ifeat].plot(plot_x, plot_y, color=cmap(norm(color_val)), ls='--') 
    axs[ifeat].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    # ============================================================================================
    # Plot
    # ============================================================================================
                
    # Plot MTNG
    for snap in MTNG_snaps:
        MTNG_z, MTNG_bin_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',  
                                                        f'med_{bin_type}', feature)

        # Duplicate z to the length of the bin data
        MTNG_z = np.repeat(MTNG_z, len(MTNG_bin_data['median']))
        # Plot
        for i in range(len(MTNG_z)):
            axs[ifeat].errorbar(MTNG_z[i], MTNG_feat['median'][i],
                            yerr=[[MTNG_feat['median'][i]-MTNG_feat['min'][i]], 
                                [MTNG_feat['max'][i]-MTNG_feat['median'][i]]],
                        color=cmap(norm(MTNG_bin_data['median'][i])), fmt='.')
    
    # Final edit
    # if feature == 'depth':
    #     Y_label = r"$\mathcal{D}$"
    #     if bin_type == 'peakHeight':
    #         axs.set_ylim(2.3, 4)
    #     elif bin_type == 'NFWconc':
    #         axs.set_ylim(2.3, 3.5)
    # elif feature == 'width_dimless':
    #     Y_label = r"$\mathcal{W}$"
    # elif feature == 'DWratio':
    #     Y_label = r"$\mathcal{D}/\mathcal{W}$"
    #     if bin_type == 'NFWconc':
    #         axs.set_ylim(1, 3.8)
    axs[ifeat].set_ylabel(Ylabels[ifeat])
    # axs.legend(loc='best')
    
# ============================================================================================
# Save the plot
# ============================================================================================
axs[1].set_xlabel('z')
save_dir = f'result/paper_plots/fig4/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig4_{bin_type}_{simu}')
plt.close()