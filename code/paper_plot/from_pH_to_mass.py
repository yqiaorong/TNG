import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
from colossus.cosmology import cosmology
from colossus.lss import peaks
plt.style.use('code/style.mplstyle')
from func import load_stats, load_all_z, plot_feature, fitting
cosmology.setCosmology('planck15')


x_type = 'mass'
print('')
print(f'>>> Plot depth and width vs {x_type} <<<')
print('')

root_dir = f'result/bootstrap_stats_DK14/with_{x_type}/'
simu = 'Hydro'
features = ['width_dimless']
Ylabels = [r"$\mathcal{W}$"]


# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================

def mass_width(inputs, a, b, c):
    mass, zval = inputs
    return a*np.log10(mass)**b / (zval + 1)**c
    # return a*np.exp(b*np.log10(mass)) / (zval + 1)**c
def width_pH_mass(inputs, a, b):
    mass, zval = inputs
    peakHeight = peaks.peakHeight(mass, zval)
    return a*peakHeight**b

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(1, 1, figsize=(4, 3.5), dpi=500, sharex=True, constrained_layout=True)


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
                ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# -----------------------------------------------------------------------------------------

# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.ax.tick_params(axis='x', rotation=0) 
cb.set_label('z')

# ============================================================================================
# Plot
# ============================================================================================
for ifeat, feature in enumerate(features):
    
    all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
    all_y_med, all_y_min, all_y_max = [], [], []
        
    # Plot MTNG
    for snap in MTNG_snaps:
        MTNG_z, MTNG_mass, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                    'med_mass', feature)
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
    popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(mass_width, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = [x_type, 'z', feature],
                                                            plot_info = [axs, cmap, norm, '--'], 
                                                            bootstrap=True)
    # print(popt, red_chi2)
    # print(perr)
    print(feature)
    print('')
    
    # Fit the data
    popt, perr, red_chi2, y_fit, axs, fitted_data = fitting(width_pH_mass, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = [x_type, 'z', feature],
                                                            plot_info = [axs, cmap, norm, 'dotted'], 
                                                            bootstrap=True)
    # print(popt, red_chi2)
    # print(perr)
    print(feature, 'from ph to mass')
    print('')
    
    # Save fitted data
    fitted_data_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
    if not os.path.exists(fitted_data_dir):
        os.makedirs(fitted_data_dir)
    np.save(fitted_data_dir + f'{x_type}_{feature}_fitted_data.npy', fitted_data)
    
    axs.set_xscale('log')
    axs.set_ylabel(Ylabels[ifeat])
    # axs[ifeat].legend(loc='best')
    
# ============================================================================================
# Save the plot
# ============================================================================================

axs.set_xlabel(r"$M_{200m} / M_\odot$")
  
save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/width_pH_mass_{simu}')
plt.close()