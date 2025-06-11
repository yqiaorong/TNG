import os
import numpy as np
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import load_stats, plot_feature, fitting

x_type = 'formzOLD'
print('')
print(f'>>> Plot depth and width vs {x_type} <<<')
print('')

root_dir = f'result/bootstrap_stats_DK14/with_{x_type}/'
simu = 'Hydro'
features = ['depth', 'width_dimless']
Ylabels = [r"$\mathcal{D}$", r"$\mathcal{W}$"]

# ============================================================================================
# Define the fitting function with two variables
# ============================================================================================

def depth(inputs, a, b):
    x, zval = inputs
    return a*x + b

# def depth_from_mass(inputs, a):
#     x, z = inputs
#     return 0.008 * (a-np.log10(x))**2.23
def depth_from_mass(inputs, a, A, B):
    x, z = inputs
    return A * (a-np.log10(x))**B

def width(inputs, a, b):
    x, zval = inputs
    return a*np.log(x) + b

# def width_from_mass(inputs, a):
#     x, z = inputs
#     return 4.9*1E6 * (a-np.log10(x))**(-5.6)
def width_from_mass(inputs, A, a, B):
    x, z = inputs
    return A * (a-np.log10(x))**B
    
# def DW(inputs, a, b, c):
#     x, zval = inputs
    # if x_type == 'formzOLD':
    #     return (-0.34200887*x+3.14556083)/(0.5975754*np.log(x)+2.04735604) # formzOLD
    # if x_type == 'formzOLD':
    #     return (-0.59763004*x+3.34975648)/(0.76350522*np.log(x)+2.10171115) # formzOLD
    
# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(2, 1, figsize=(4, 6), dpi=500, sharex=True, constrained_layout=True)

# ============================================================================================
# Set up the colorbar
# ============================================================================================

MTNG_snaps = [264]
MTNG_dir = f'{root_dir}/MTNG/{simu}-Arepo/MTNG-L500-4320-A/Nboots_1024/'

# all_z = load_all_z(MTNG_dir, MTNG_snaps)
# all_z.append(2.1)
all_z = np.array([0, 1, 2])

# Set up the colorbar
# -----------------------------------------------------------------------------------------
cmap = plt.get_cmap('viridis', len(all_z))
norm = BoundaryNorm(all_z, cmap.N)
# cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
#                 ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# # -----------------------------------------------------------------------------------------

# ============================================================================================
# Plot
# ============================================================================================
for ifeat, feature in enumerate(features):
    print(feature)
    if feature == 'depth':
        fit_func = depth
    elif feature == 'width_dimless':
        fit_func = width
    
    all_z, all_x_med, all_x_min, all_x_max = [], [], [], []
    all_y_med, all_y_min, all_y_max = [], [], []
        
    # Plot MTNG
    for snap in MTNG_snaps:
        MTNG_z, MTNG_x_data, MTNG_feat = load_stats(MTNG_dir, f'snap_{snap}_Rsp_stats.npy',
                                                    f'med_{x_type}', feature)
        
        # select x data which >= 0.25
        MTNG_feat = {k: v[MTNG_x_data['median'] >= 0.25] for k, v in MTNG_feat.items()}
        MTNG_x_data = {k: v[MTNG_x_data['median'] >= 0.25] for k, v in MTNG_x_data.items()}
        
        plot_feature(simu, np.round(MTNG_z, 1), MTNG_x_data, MTNG_feat, [axs[ifeat], cmap, norm])
        
        # Append the data
        all_x_med.append(MTNG_x_data['median'])
        all_x_min.append(MTNG_x_data['min'])
        all_x_max.append(MTNG_x_data['max'])
        all_y_med.append(MTNG_feat['median'])
        all_y_min.append(MTNG_feat['min'])
        all_y_max.append(MTNG_feat['max'])
        # Duplicate z to the same length as the data
        all_z.append([MTNG_z]*len(MTNG_x_data['median']))
    
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
        
        # Save fitted data
        fitted_data_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
        if not os.path.exists(fitted_data_dir):
            os.makedirs(fitted_data_dir)
        np.save(fitted_data_dir + f'{x_type}_{feature}_fitted_data.npy', fitted_data)
        
    # ============================================================================================
    # Plot the fitting from mass equation
    # ============================================================================================
    if feature == 'depth':
        popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(depth_from_mass, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = [x_type, 'z', feature],
                                                            plot_info = [axs[ifeat], cmap, norm,  'dotted'], 
                                                            bootstrap=False)
    if feature == 'width_dimless':
        popt, perr, red_chi2, y_fit, axs[ifeat], fitted_data = fitting(width_from_mass, 
                                                            values = [all_x_med, all_z, all_y_med, all_y_min, all_y_max],
                                                            labels = [x_type, 'z', feature],
                                                            plot_info = [axs[ifeat], cmap, norm,  'dotted'], 
                                                            bootstrap=False, 
                                                            p0 = [5*10**6, 1, -1],
                                                            bounds=([4*10**6, 0, -10], [5.5*10**6, 10, 0])
                                                            )
    print(popt)
    print(perr)           
    print('')
    
    axs[ifeat].set_ylabel(Ylabels[ifeat])
    # axs.legend(loc='best')
    
# ============================================================================================
# Save the plot
# ============================================================================================
axs[-1].set_xlabel(r'$z_{\rm form}$')
save_dir = f'result/paper_plots/fig3_fit/MTNG-Hydro/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig3_{x_type}_{simu}')
plt.close()