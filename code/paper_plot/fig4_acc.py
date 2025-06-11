import os
import numpy as np
from func import get_bins, load_stats, fitting
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')

bin_type = 'accretions'
label = r'$\Gamma$'
min_bin, max_bin, bin_width = 0, 6, 0.5
num_bins, all_bins, _, _ = get_bins(min_bin, max_bin, bin_width)
  
print('')
print(f'>>> Plot depth and width vs z ({bin_type}) <<<')
print('')

def accret_depth0(inputs, a, b, c):
    zval, accret = inputs    
    return a*zval**b + c

def accret_depth(inputs, a, b, c, d):
    zval, accret = inputs    
    return a*zval**b + c*accret**d
    
def accret_depth2(inputs, a, b):
    zval, accret = inputs 
    return a*(zval+1)**b 

def accret_width0(inputs, a, b, c):
    zval, accret = inputs   
    return a*zval**b + c

def accret_width(inputs, a, b, c, d):
    zval, accret = inputs   
    return a*zval**b + c*accret**d

def accret_width2(inputs, a, b):
    zval, accret = inputs 
    return a*(zval+1)**b 
    
# def accret_DW(inputs, a, b, c):
#     zval, accret = inputs   
#     return a*zval**b + c


root_dir = f'result/bootstrap_stats/with_{bin_type}/'
simu = 'Hydro'
features = ['depth', 'width_dimless',]
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$",]

fig, axs = plt.subplots(2, 1, figsize=(4, 6), dpi=500, sharex=True, constrained_layout=True) 


# Set up the colorbar
cmap = plt.get_cmap('plasma', len(all_bins))
if bin_type == 'mass':
    bound = np.logspace(min(all_bins), max(all_bins)+0.01, len(all_bins)+1) 
else:
    bound = np.linspace(min(all_bins), max(all_bins)+0.01, len(all_bins)+1) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                    ax=axs[-1], orientation='horizontal', spacing='proportional', ticks=bound)

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.ax.tick_params(axis='x', rotation=0) 

cb.set_label(label)

# ============================================================================================
# Load data 
# ============================================================================================

MTNG_snaps = [264, 237, 214, 179, 151, 129]
MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'      

fit_func0 = [accret_depth0, accret_width0]
fit_func1 = [accret_depth, accret_width]
fit_func2 = [accret_depth2, accret_width2]

for ifeat, feature in enumerate(features):

    print(feature)

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
            axs[ifeat].errorbar(MTNG_z[i], MTNG_feat['median'][i],
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
    popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(fit_func0[ifeat], 
                                                            values = [all_z, all_x_med, all_y_med, all_y_min, all_y_max],
                                                            labels = ['z', bin_type, feature],
                                                            plot_info = [axs[ifeat], None, norm, '--'], 
                                                            bootstrap=True)
    print(popt, red_chi2)
    print(perr)
    popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(fit_func1[ifeat], 
                                                            values = [all_z, all_x_med, all_y_med, all_y_min, all_y_max],
                                                            labels = ['z', bin_type, feature],
                                                            plot_info = [axs[ifeat], cmap, norm, 'dotted'], 
                                                            bootstrap=True)
    print(popt, red_chi2)
    print(perr)
    popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(fit_func2[ifeat], 
                                                            values = [all_z, all_x_med, all_y_med, all_y_min, all_y_max],
                                                            labels = ['z', bin_type, feature],
                                                            plot_info = [axs[ifeat], None, norm, 'solid'], 
                                                            bootstrap=True)
    print(popt, red_chi2)
    print(perr)
    print('')

    axs[ifeat].set_ylabel(Ylabels[ifeat])
    # axs[ifeat].legend(loc='best')
    
# ============================================================================================
# Save the plot
# ============================================================================================

axs[-1].set_xlabel('z')
save_dir = f'result/paper_plots/fig4_fit/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig4_{bin_type}_{simu}')
plt.close()