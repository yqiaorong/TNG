import numpy as np
from matplotlib import pyplot as plt

def load_stats_key(dir, fname, key):
    try:
        data = np.load(f'{dir}/{fname}', allow_pickle=True).item()
    except ValueError: 
        data = np.load(f'{dir}/{fname}', allow_pickle=True)
        data = data[0]
    
    return data[key]

    
def load_stats(dir, fname, xlabel, ylabel=None, bin_val='z', Pick_z=True):
    try:
        data = np.load(f'{dir}/{fname}', allow_pickle=True).item()
    except ValueError: 
        data = np.load(f'{dir}/{fname}', allow_pickle=True)
        data = data[0]
    
    if Pick_z == True:
        z = np.round(data[bin_val], 3)
    else:
        z = 0

    x, xmin, xmax = data[xlabel][1], data[xlabel][0], data[xlabel][2]
    if xlabel == 'med_mass':
        x = np.array([10**10*m for m in x])
        xmin = np.array([10**10*m for m in xmin])
        xmax = np.array([10**10*m for m in xmax])
    x_dict = {'median': x, f'min': xmin, f'max': xmax}
     
    if ylabel is not None:
        y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]
        y_dict = {'median': y, f'min': ymin, f'max': ymax}
    else:
        y_dict = None
    
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


def plot_theory_feature(z, x_dict, y_dict, plot_info, label=None, ls='--'):
    axs, cmap, norm = plot_info
                     
    if label==None:
        axs.plot(x_dict, y_dict, color=cmap(norm(z)), ls=ls)
    else:
        axs.plot(x_dict, y_dict, color=cmap(norm(z)), ls=ls, label=label)
         
               
def plot_slope_profile(radius, slopes, save_dir, save_name):
    fig_s, axs_s = plt.subplots(1, 1, figsize=(4, 3.3), dpi=500, sharex=True, constrained_layout=True)
    axs_s.plot(radius, slopes)
    axs_s.set_xlabel(r'$r/R_{200}$')
    axs_s.set_ylabel(r'$\nabla \rho/\rho$')
    axs_s.set_xscale('log')
    axs_s.set_ylim(-6, 0)
    plt.savefig(f'{save_dir}/{save_name}.png')
    plt.close(fig_s)
        
        
# Define the fitting with two variables
def fitting(func, values, labels, plot_info, p0=None, bounds=None, bootstrap=False):
    from scipy.optimize import curve_fit
    
    x1, x2, y, ymax, ymin = np.asarray(values[0]), np.asarray(values[1]), np.asarray(values[2]),\
        np.asarray(values[3]), np.asarray(values[4])
    xscales, bins, features = labels[0], labels[1], labels[2]

    # Get the mean eroors
    yerr_mean = (ymax - ymin) / 2

    # yerr plus epsilon if yerr is zero
    yerr_mean[yerr_mean == 0] = 1e-1
    
    # Curve fitting with input sigma
    if p0 is None:
        p0 = [1] * (func.__code__.co_argcount - 1) 
    # For fiiting width as a function of z
    # p0= [-7.20882095, 1.63511888, -7.05248173E-2, -1.17051206E1,
    #      -3.66871643E-1,  1.07646616,  1.09160168, -9.71760263e-05]
    
    if bounds is None:
        popt, pcov = curve_fit(func, (x1, x2), y, p0=p0,   
                           sigma=yerr_mean, 
                           maxfev=1000000)
    else:
        popt, pcov = curve_fit(func, (x1, x2), y, p0=p0,   
                            sigma=yerr_mean, 
                            #    bounds = ([-5] * len(p0), [5] * len(p0)),
                            bounds = bounds,
                            maxfev=1000000)
    print('Optimal parameters:', popt)
    if bootstrap:
        perr = boots_err([x1, x2], y, yerr_mean, func, popt)
    else:
        perr = np.sqrt(np.diag(pcov))
        print('Errors:', perr)

    # Compute reduced chi-square 
    y_fit = func([x1, x2], *popt)
    if np.all(popt == 1):
        red_chi2 = np.sum((y-y_fit)**2 / yerr_mean**2) / len(y)
    else:
        red_chi2 = np.sum((y-y_fit)**2 / yerr_mean**2) / (len(y) - len(popt))
    print('Non reduced chi2:', np.sum((y-y_fit)**2 / yerr_mean**2))
    print('Reduced chi2:', red_chi2)
    
    # Plot the fitting
    axs, cmap, norm, linestyle = plot_info      

    # fitted data
    fitx1, fitx2, fity = [], [], []
    for uniq_x2 in np.unique(x2):
        if cmap is None:
            color= 'black'
        else:
            color = cmap(norm(uniq_x2))
        ib = np.where(uniq_x2 == x2)[0]
        
        # Find minimum and maximum values
        # ---------------------------------------
        x1_fit = np.linspace(0, 2, 10)
        # x1_fit = np.logspace(13, 15.5, 6)
        x2_fit = np.repeat(uniq_x2, len(x1_fit))
        # ---------------------------------------
        
        y_fit = func([x1_fit, x2_fit], *popt)
        if 0 in ib:
            axs.plot(x1_fit, y_fit, # c='black', 
                     c=color, 
                     ls=linestyle, label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}')
        else:
            axs.plot(x1_fit, y_fit, # c='black', 
                     c=color, 
                     ls=linestyle)
        fitx1.append(x1_fit)
        fitx2.append(np.repeat(uniq_x2, len(x1_fit)))
        fity.append(y_fit)

    fitted_data = {
                    xscales:    np.asarray(fitx1),
                    bins:       np.asarray(fitx2),
                    features:   np.asarray(fity),
                    'red_chi2': red_chi2
                  }
            
    return popt, perr, red_chi2, y_fit, axs, fitted_data


def boots_err(xdata, ydata, yerr, your_model, popt):
    from scipy.optimize import curve_fit
    
    x1, x2 = xdata[0], xdata[1]
    
    # Bootstrap
    n_bootstrap = 100
    params_bootstrap = []

    for _ in range(n_bootstrap):
        indices = np.random.choice(len(ydata), len(ydata), replace=True)
        x1_bs, x2_bs = x1[indices], x2[indices]
        y_bs, yerr_bs = ydata[indices], yerr[indices]
        
        try:
            popt_bs, _ = curve_fit(your_model, (x1_bs, x2_bs), y_bs, sigma=yerr_bs, absolute_sigma=True)
            params_bootstrap.append(popt_bs)
        except RuntimeError:
            continue  # skip failed fits

    params_bootstrap = np.array(params_bootstrap)
    lower_bounds = np.percentile(params_bootstrap, 16, axis=0)
    upper_bounds = np.percentile(params_bootstrap, 84, axis=0)

    # asymmetric errors
    errors_minus = popt - lower_bounds
    errors_plus = upper_bounds - popt
    # print('errors_minus:', errors_minus)
    # print('errors_plus:', errors_plus)
    # Get average error
    errors_avg = (errors_minus + errors_plus) / 2
    print('average error:', errors_avg)
    return [errors_plus, errors_minus]

    
def density_gradient_profile(
    r: np.ndarray,
    rho_s: float,
    r_s: float,
    alpha: float,
    r_t: float,
    beta: float,
    gamma: float,
    rho_g: float,
    b_e: float,
    S_e: float,
    R_200_mean: float,) -> np.ndarray:

    rho_inner = rho_s * np.exp(-(2.0 / alpha) * (np.power(r / r_s, alpha) - 1.0))
    f_trans = np.power(1.0 + np.power(r / r_t, beta), -gamma / beta)
    rho_outer = rho_g * (b_e * np.power(r / (5.0 * R_200_mean), -S_e) + 1.0)

    rho = rho_inner * f_trans + rho_outer

    d_rho_inner_dr = - (2.0 / r_s) * np.power(r / r_s, alpha - 1) * rho_inner
    d_ftrans_dr = np.power(
        1.0 + np.power(r / r_t, beta), -(1.0 + gamma / beta)
    ) * (- gamma / r_t) * np.power(r / r_t, beta - 1.0)

    d_rho_outer_dr = - (rho_g * b_e * S_e) / (5.0 * R_200_mean) * np.power(r / (5.0 * R_200_mean), -(S_e + 1.0))

    d_rho_dr = d_rho_inner_dr * f_trans + rho_inner * d_ftrans_dr + d_rho_outer_dr

    return (r / rho) * d_rho_dr
    
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
        try:
            data = np.load(f'{Dir}/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        except ValueError: 
            data = np.load(Dir+f'snap_{snap}_Rsp_stats.npy', allow_pickle=True)
            data = data[0]
        all_z.append(np.round(data['z'], 2))
    return all_z



# def delta_c(z):
#     """Compute the characteristic overdensity at redshift z."""
#     import astropy.units as u
#     from astropy.constants import M_sun, G
#     from astropy.cosmology import FlatLambdaCDM
    
#     # Define own cosmology
#     h = 0.6774
#     Om0 = 0.3089
#     Ob0 = 0.0486
#     cosmo = FlatLambdaCDM(H0 = h * 100 * u.km / u.s / u.Mpc, 
#                             Om0=Om0, Ob0=Ob0, Tcmb0=2.725)
    
#     # Compute the critical density
#     H = cosmo.H(z) # km / (Mpc s)
#     H = H.to(u.m / u.s / u.m)
#     rho_c = 3 * H**2 / ( 8 * np.pi * G ) 
#     rho_c = rho_c.to(M_sun / u.kpc**3) # Msun / kpc**3
    
#     # Compute the characteristic overdensity
#     delta_rho = rho_c * 1.686
    
#     return delta_rho

# def character_mass(char_r, char_rho):
    return 4/3 * np.pi * char_r**3 * char_rho