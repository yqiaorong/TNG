import numpy as np
from matplotlib import pyplot as plt
    
    
def load_stats(dir, fname, xlabel, ylabel, bin_val='z'):
    data = np.load(f'{dir}/{fname}', allow_pickle=True).item()

    z = np.round(data[bin_val], 3)
    x, xmin, xmax = data[xlabel][1], data[xlabel][0], data[xlabel][2]
    print(x)
    if xlabel == 'med_mass':
        x = np.array([10**10*m for m in x])
        xmin = np.array([10**10*m for m in xmin])
        xmax = np.array([10**10*m for m in xmax])
    y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]

    x_dict = {'median': x, f'min': xmin, f'max': xmax}
    y_dict = {'median': y, f'min': ymin, f'max': ymax}
    return z, x_dict, y_dict

def plot_feature(simu, z, x_dict, y_dict, plot_info, label=None, ls='.'):
    axs, cmap, norm = plot_info
                     
    if label==None:
        axs.errorbar(x_dict['median'], y_dict['median'],
                    xerr=[x_dict['median']-x_dict['min'], x_dict['max']-x_dict['median']],
                    yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
                    color=cmap(norm(z)), fmt=ls)
    else:
        axs.errorbar(x_dict['median'], y_dict['median'],
                    xerr=[x_dict['median']-x_dict['min'], x_dict['max']-x_dict['median']],
                    yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
                    color=cmap(norm(z)), fmt=ls, label=label)
        
        
# Define the fitting with two variables
def fitting(bin_type, feature, bins, x, y, ymin, ymax, func, plot_info, xmin=None, xmax=None):
    from scipy.odr import Model, RealData, ODR
    from scipy.optimize import curve_fit
    
    bins = np.asarray(bins)

    # Get the mean eroors
    if xmin is not None:
        x_err = (xmax - xmin) / 2
    y_err = (ymax - ymin) / 2

    # yerr plus epsilon if yerr is zero
    y_err[y_err == 0] = 1e-1
    
    # Curve fitting with input sigma
    p0 = [1] * (func.__code__.co_argcount - 1) 
    # For fiiting width as a function of z
    # p0= [-7.20882095, 1.63511888, -7.05248173E-2, -1.17051206E1,
    #      -3.66871643E-1,  1.07646616,  1.09160168, -9.71760263e-05]

    popt, pcov = curve_fit(func, (x, bins), y, p0=p0,   
                           sigma=y_err, 
                           maxfev=1000000)
    perr = np.sqrt(np.diag(pcov))

    # Compute reduced chi-square 
    y_fit = func([x, bins], *popt)
    red_chi2 = np.sum((y-y_fit)**2 / y_err**2) / (len(y) - len(popt))
    
    # Plot the fitting
    axs, cmap, norm = plot_info
    
    # fitted data
    fitx, fity, fitbin = [], [], []
    for uniq_bin in np.unique(bins):
        ib = np.where(uniq_bin == bins)[0]
        if 0 in ib:
            # sort according to x
            idx = np.argsort(x[ib])
            axs.plot(x[ib][idx], y_fit[ib][idx], c=cmap(norm(uniq_bin)), label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}', ls='--')
        else:
            idx = np.argsort(x[ib])
            axs.plot(x[ib][idx], y_fit[ib][idx], c=cmap(norm(uniq_bin)), ls='--')
        fitx.append(x[ib][idx])
        fity.append(y_fit[ib][idx])
        fitbin.append(np.repeat(uniq_bin, len(x[ib][idx])))

    fitted_data = {
        bin_type: np.concatenate(fitx),
        feature: np.concatenate(fity),
        'z': np.concatenate(fitbin),
        'red_chi2': red_chi2
    }
   
    return popt, perr, red_chi2, y_fit, axs, fitted_data
    
    
# def load_stats_per_bin(dir, fname_format, fname_val, 
#                        xlabel, ylabel, bin_name, bin_start_val):
    
#     all_x, all_y, all_y_min, all_y_max = [], [], [], []
    
#     for val in fname_val:
#         data = np.load(dir+fname_format.format(val), allow_pickle=True).item()

#         bins = data[bin_name]
#         bin_start_idx =  np.where(bins == bin_start_val)[0]
        
#         if len(bin_start_idx) == 0:
#             pass
#         else:
#             x = np.round(data[xlabel], 3)
            
#             y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]
            
#             # Select the corresponding bin vals
#             y, ymin, ymax = y[bin_start_idx], ymin[bin_start_idx], ymax[bin_start_idx]
           
#             # Append data
#             all_y.append(y)
#             all_y_min.append(ymin)
#             all_y_max.append(ymax)
#             all_x.append(x)
    
#     if len(all_x) < 2:
#         all_y = np.array(all_y).ravel()
#         all_y_min = np.array(all_y_min).ravel()
#         all_y_max = np.array(all_y_max).ravel()
#     else:
#         all_y = np.concatenate(all_y)
#         all_y_min = np.concatenate(all_y_min)
#         all_y_max = np.concatenate(all_y_max)

#     # Convert the data to dict
#     y_dict = {'median': all_y, f'min': all_y_min, f'max': all_y_max}

#     return all_x, y_dict

# def plot_feature_vs_z(simu, bin_val, all_x, y_dict, plot_info):
   
#     if len(all_x) == 0:
#         pass
#     else:
#         axs, cmap, norm = plot_info
#         if simu == 'DM':
#             ls = '--'
#         else:
#             ls = '-'
#         axs.errorbar(all_x, y_dict['median'],
#                      yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
#                      color=cmap(norm(bin_val)), fmt='.')

def get_bins(min_bin, max_bin, bin_width):
    num_bins = int((max_bin - min_bin)/bin_width)
    all_bins = np.arange(min_bin, max_bin+bin_width, bin_width)
    bins_starts = all_bins[:-1]
    bins_ends = all_bins[1:]
    return num_bins, all_bins, bins_starts, bins_ends

def load_all_z(Dir, snaps):
    all_z = []
    for snap in snaps:
        data = np.load(Dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        all_z.append(np.round(data['z'], 2))
    return all_z

def delta_c(z):
    """Compute the characteristic overdensity at redshift z."""
    import astropy.units as u
    from astropy.constants import M_sun, G
    from astropy.cosmology import FlatLambdaCDM
    
    # Define own cosmology
    h = 0.6774
    Om0 = 0.3089
    Ob0 = 0.0486
    cosmo = FlatLambdaCDM(H0 = h * 100 * u.km / u.s / u.Mpc, 
                            Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    
    # Compute the critical density
    H = cosmo.H(z) # km / (Mpc s)
    H = H.to(u.m / u.s / u.m)
    rho_c = 3 * H**2 / ( 8 * np.pi * G ) 
    rho_c = rho_c.to(M_sun / u.kpc**3) # Msun / kpc**3
    
    # Compute the characteristic overdensity
    delta_rho = rho_c * 1.686
    
    return delta_rho

def character_mass(char_r, char_rho):
    return 4/3 * np.pi * char_r**3 * char_rho