import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, load_all_z, plot_feature, fitting

x_type = 'peakHeight'
print('')
print(f'>>> Plot depth and width vs {x_type} <<<')
print('')

root_dir = f'result/bootstrap_stats_DK14/with_{x_type}/'
simu = 'Hydro'
features = ['depth', 'width_dimless', 'DWratio']
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$", r"$\mathcal{D}/\mathcal{W}$"]

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================

def depth(inputs, a, b, c, d):
    x, zval = inputs
    return a*x**b + c*zval**d

def width(inputs, a, b, c):
    x, zval = inputs
    return a*x**b + c*zval

def DW(inputs, a, b, c, d):
    x, zval = inputs
    return a*x**b + c*zval**d

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(3, 1, figsize=(4, 9), dpi=500, sharex=True, constrained_layout=True)


MTNG_snaps = [264, 237, 214, 179, 151, 129]
MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

all_z = load_all_z(MTNG_dir, MTNG_snaps)
all_z.append(2.1)

# ============================================================================================
# Set up the colorbar
# ============================================================================================
# -----------------------------------------------------------------------------------------
cmap = plt.get_cmap('viridis', len(MTNG_snaps))
bound = all_z
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                ax=axs[2], orientation='horizontal', spacing='proportional', ticks=bound)
# -----------------------------------------------------------------------------------------

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.ax.tick_params(axis='x', # rotation=70
                  ) 
cb.set_label('z')

# ============================================================================================
# Plot
# ============================================================================================
for ifeat, feature in enumerate(features):
    print(feature)
    if feature == 'depth':
        fit_func = depth
    elif feature == 'width_dimless':
        fit_func = width
    elif feature == 'DWratio':
        fit_func = DW      

    all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
    all_y_med, all_y_min, all_y_max = [], [], []
            
    # Plot MTNG
    for snap in MTNG_snaps:
        MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                    f'med_{x_type}', feature)
        plot_feature(simu, np.round(MTNG_z, 1), MTNG_mass, MTNG_feat, [axs[ifeat], cmap, norm])
        
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
    if fit_func is not None:
        popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(fit_func, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = [x_type, 'z', feature],
                                                            plot_info = [axs[ifeat], cmap, norm, '--'], 
                                                            bootstrap=True)
        print(popt)
        print(perr)
        print(red_chi2)
        print('')
        
        # Save fitted data
        fitted_data_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
        if not os.path.exists(fitted_data_dir):
            os.makedirs(fitted_data_dir)
        np.save(fitted_data_dir + f'{x_type}_{feature}_fitted_data.npy', fitted_data)
    
    axs[ifeat].set_ylabel(Ylabels[ifeat])
    # axs[ifeat].legend(loc='best')
    
# ============================================================================================
# Save the plot
# ============================================================================================

axs[0].set_ylim(2.25, 4)
axs[-1].set_xlabel(r"$\nu$")

save_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig3_{x_type}_{simu}')
plt.close()