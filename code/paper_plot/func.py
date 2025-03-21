import numpy as np
from matplotlib import pyplot as plt
    
def load_stats(dir, fname, xlabel, ylabel, bin_val='z'):
    data = np.load(f'{dir}/{fname}', allow_pickle=True).item()

    z = np.round(data[bin_val], 3)
    x, xmin, xmax = data[xlabel][1], data[xlabel][0], data[xlabel][2]
    if xlabel == 'med_mass':
        x = np.array([10**10*m for m in x])
        xmin = np.array([10**10*m for m in xmin])
        xmax = np.array([10**10*m for m in xmax])
    y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]

    x_dict = {'median': x, f'min': xmin, f'max': xmax}
    y_dict = {'median': y, f'min': ymin, f'max': ymax}
    return z, x_dict, y_dict

# Define the fitting with two variables
def fitting(bins, x, y, ymin, ymax, func, plot_info, xmin=None, xmax=None):
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

    for uniq_bin in np.unique(bins):
        ib = np.where(uniq_bin == bins)[0]
        if 0 in ib:
            # sort according to x
            idx = np.argsort(x[ib])
            axs.plot(x[ib][idx], y_fit[ib][idx], c=cmap(norm(uniq_bin)), label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}', ls='--')
        else:
            idx = np.argsort(x[ib])
            axs.plot(x[ib][idx], y_fit[ib][idx], c=cmap(norm(uniq_bin)), ls='--')
   
    return popt, red_chi2, y_fit, axs
    
def load_stats_per_bin(dir, fname_format, fname_val, 
                       xlabel, ylabel, bin_name, bin_start_val):
    
    all_x, all_y, all_y_min, all_y_max = [], [], [], []
    
    for val in fname_val:
        data = np.load(dir+fname_format.format(val), allow_pickle=True).item()

        bins = data[bin_name]
        bin_start_idx =  np.where(bins == bin_start_val)[0]
        
        if len(bin_start_idx) == 0:
            pass
        else:
            x = np.round(data[xlabel], 3)
            
            y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]
            
            # Select the corresponding bin vals
            y, ymin, ymax = y[bin_start_idx], ymin[bin_start_idx], ymax[bin_start_idx]
           
            # Append data
            all_y.append(y)
            all_y_min.append(ymin)
            all_y_max.append(ymax)
            all_x.append(x)
    
    if len(all_x) < 2:
        all_y = np.array(all_y).ravel()
        all_y_min = np.array(all_y_min).ravel()
        all_y_max = np.array(all_y_max).ravel()
    else:
        all_y = np.concatenate(all_y)
        all_y_min = np.concatenate(all_y_min)
        all_y_max = np.concatenate(all_y_max)

    # Convert the data to dict
    y_dict = {'median': all_y, f'min': all_y_min, f'max': all_y_max}

    return all_x, y_dict

def plot_feature(simu, z, x_dict, y_dict, plot_info):
    axs, cmap, norm = plot_info
    if simu == 'DM':
        ls = '--'
    else:
        ls = '-'
                    
    # axs.plot(x_dict['median'], y_dict['median'], color=cmap(norm(z)), lw=1, alpha=0.5, linestyle=ls) 
    axs.errorbar(x_dict['median'], y_dict['median'],
                    xerr=[x_dict['median']-x_dict['min'], x_dict['max']-x_dict['median']],
                    yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
                    color=cmap(norm(z)), fmt='.')
        
def plot_feature_vs_z(simu, bin_val, all_x, y_dict, plot_info):
   
    if len(all_x) == 0:
        pass
    else:
        axs, cmap, norm = plot_info
        if simu == 'DM':
            ls = '--'
        else:
            ls = '-'
        axs.errorbar(all_x, y_dict['median'],
                     yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
                     color=cmap(norm(bin_val)), fmt='.')

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

# Fitting

def x_z_fit_d(X, redshifts, Y, Yerr, axs, cmap, norm):
    from scipy.optimize import curve_fit
    
    def func(Inputs, a, b, c, d, e, f, g, h, i):
        """Params:
            Inputs: (x, redshifts)
        """
        x, z = Inputs
        return a*x + b*z + d + c*x/(z-e) + f*z*(x-g*z)**2 + h*z/(x-i*z)

    popt, pcov = curve_fit(func, (X, redshifts), Y, p0=[1]*9, maxfev=10000)
    
    Y_fit = func((X, redshifts), *popt)

    # Calculate the reduced chi-square
    if np.ndim(Yerr) == 2:
       Yerr = np.mean(Yerr, axis=0)
       
    # Compute effective errors
    red_chi2 = np.sum((Y - Y_fit)**2 / Yerr**2) / (len(Y) - len(popt))
    # Plot
    for z in np.unique(redshifts):
        iz = np.where(redshifts == z)[0]
        if 0 in iz:
            axs.plot(X[iz], Y_fit[iz], c=cmap(norm(np.round(z, 2))), label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}', ls='--')
        else:
            axs.plot(X[iz], Y_fit[iz], c=cmap(norm(np.round(z, 2))), ls='--')

    # Confidence interval
    perr = np.sqrt(np.diag(pcov))

    dof = max(1, len(Y)-len(popt)) 
    
    # 99% confidence level
    from scipy.stats import t
    alpha = 0.01 
    t_score = t.ppf(1 - alpha/2, dof)

    # Confidence interval = popt ± (t-score * std_err)
    ci_lower = popt - t_score * perr
    ci_upper = popt + t_score * perr

    # Print results
    for i, (p, lo, up) in enumerate(zip(popt, ci_lower, ci_upper)):
        print(f"Parameter {i}: {p:.4f} (99% CI: {lo:.4f} to {up:.4f})")

def x_z_fit_w(X, redshifts, Y,Yerr, axs, cmap, norm):
    from scipy.optimize import curve_fit
    
    def func(Inputs, a, b, d, A, D, E, F):
        """Params:
            Inputs: (x, redshifts)
        """
        x, z = Inputs
        return a*x + b*z + d + A*x**2 + D*z**3 + F*x**2*z + E*x*z**2 # + G*z/x  
    
    # def func(Inputs, a, d, A, G):
    #     """Params:
    #         Inputs: (x, redshifts)
    #     """
    #     x, z = Inputs
    #     return a*x + d + A*x**2 + G/x 
    
    # def func(Inputs, a, b, c, d, e):
    #     """Params:
    #             Inputs: (x, redshifts)
    #     """
    #     x, z = Inputs
    #     return a*z/x*np.tanh(-b*z*(x-d)) + c*np.tanh(-e*x)
    
    # def func(Inputs, a, b, c, e):
    #     """Params:
    #         Inputs: (x, redshifts)
    #     """
    #     x, z = Inputs
    #     return (z/x+c)*np.exp(-(a*(x-e)**2)/b)+1/x

    popt, pcov = curve_fit(func, (X, redshifts), Y, p0=[1]*7, maxfev=10000)
    
    Y_fit = func((X, redshifts), *popt)

    # Calculate the reduced chi-square
    red_chi2 = np.sum((Y - Y_fit)**2 / Yerr**2) / (len(Y) - len(popt))
    # Plot
    for z in np.unique(redshifts):
        iz = np.where(redshifts == z)[0]
        if 0 in iz:
            axs.plot(X[iz], Y_fit[iz], c=cmap(norm(np.round(z, 2))), label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}', ls='--')
        else:
            axs.plot(X[iz], Y_fit[iz], c=cmap(norm(np.round(z, 2))), ls='--')
    # axs.scatter(X, Y_fit, c=cmap(norm(np.round(z, 2))), marker='x', s=20,
    #             label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}')
    
    # 99% confidence level
    from scipy.stats import t
    perr = np.sqrt(np.diag(pcov))
    dof = max(1, len(Y)-len(popt))
    
    alpha = 0.01 
    t_score = t.ppf(1 - alpha/2, dof)

    # Confidence interval = popt ± (t-score * std_err)
    ci_lower = popt - t_score * perr
    ci_upper = popt + t_score * perr

    # Print results
    for i, (p, lo, up) in enumerate(zip(popt, ci_lower, ci_upper)):
        print(f"Parameter {i}: {p:.4f} (99% CI: {lo:.4f} to {up:.4f})")
