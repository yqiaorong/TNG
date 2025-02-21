import numpy as np

def load_stats(dir, snap, xlabel, ylabel):
    data = np.load(f'{dir}/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()

    z = np.round(data['z'], 3)
    x, xmin, xmax = data[xlabel][1], data[xlabel][0], data[xlabel][2]
    y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]

    x_dict = {'median': x, f'min': xmin, f'max': xmax}
    y_dict = {'median': y, f'min': ymin, f'max': ymax}
    return z, x_dict, y_dict
    
def load_stats_per_bin(dir, snaps, bin_name, bin_start_val, ylabel):
    
    # all_x, all_x_min, all_x_max = [], [], []
    all_y, all_y_min, all_y_max = [], [], []
    all_z = []
    
    for snap in snaps:
        data = np.load(f'{dir}/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
       
        bins = data[bin_name]
        bin_start_idx =  np.where(bins == bin_start_val)[0]
        
        if len(bin_start_idx) == 0:
            pass
        else:
            z = np.round(data['z'], 3)
            
            # x, xmin, xmax = data[xlabel][1], data[xlabel][0], data[xlabel][2]
            y, ymin, ymax = data[ylabel][1], data[ylabel][0], data[ylabel][2]
            
            # Select the corresponding bin vals
            # x, xmin, xmax = x[bin_start_idx], xmin[bin_start_idx], xmax[bin_start_idx]
            y, ymin, ymax = y[bin_start_idx], ymin[bin_start_idx], ymax[bin_start_idx]
           
            # Append data
            # all_x.append(x)
            # all_x_min.append(xmin)
            # all_x_max.append(xmax)
            all_y.append(y)
            all_y_min.append(ymin)
            all_y_max.append(ymax)
            all_z.append(z)
    
    if len(all_z) < 2:
        all_y = np.array(all_y).ravel()
        all_y_min = np.array(all_y_min).ravel()
        all_y_max = np.array(all_y_max).ravel()
    else:
        all_y = np.concatenate(all_y)
        all_y_min = np.concatenate(all_y_min)
        all_y_max = np.concatenate(all_y_max)

    # Convert the data to dict
    # x_dict = {'median': np.array(all_x), f'min': np.array(all_x_min), f'max': np.array(all_x_max)}
    y_dict = {'median': all_y, f'min': all_y_min, f'max': all_y_max}

    return all_z, y_dict

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
        
def plot_feature_vs_z(simu, bin_val, all_z, y_dict, plot_info):
   
    if len(all_z) == 0:
        pass
    else:
        axs, cmap, norm = plot_info
        if simu == 'DM':
            ls = '--'
        else:
            ls = '-'
        axs.errorbar(all_z, y_dict['median'],
                        yerr=[y_dict['median']-y_dict['min'], y_dict['max']-y_dict['median']],
                        color=cmap(norm(bin_val)), fmt='.')

def load_data(dir, snaps):
    
    for isnap, snap in enumerate(snaps): # from low z to high z (present)
            
        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = data['z']
        
        mass_bins = data['mass_bins']
        mass_bins = [10**(10+m) for m in mass_bins]
        num_bins = len(mass_bins)
        
        data = data['final_results']

        # Concatenate data
        if isnap == 0:
            tot_z = [np.round(z, 3)]*num_bins
            tot_mass_cuts = mass_bins
            tot_data = data
        else:
            tot_z += [np.round(z, 3)]*num_bins
            tot_mass_cuts += mass_bins
            tot_data = np.concatenate((tot_data, data), axis=0)
    
    return np.array(tot_z), np.array(tot_mass_cuts), tot_data

# def load_z(dir, snaps):
#     data = np.load(dir+f'/snap_{min(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
#     z_i = np.round(data['z'], 2)
#     data = np.load(dir+f'/snap_{max(snaps)}_Rsp_stats.npy', allow_pickle=True).item()
#     z_f = np.round(data['z'], 2)
#     return z_i, z_f

def load_all_z(Dir, snaps):
    all_z = []
    for snap in snaps:
        data = np.load(Dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        all_z.append(np.round(data['z'], 2))
    return all_z


### Used in plot depth vs accretion rate

def load_accret(dir, width=None):
    import numpy as np
    """
    return:
        mass_cuts:  (N,)
        z:          (M,)
        accret_med: (M, N,)
        accret_std: (M, N,)"""
    
    acc = np.load(dir, allow_pickle=True).item()
    mass_cuts = acc['mass_cuts'][:-1]
    accret_med = acc['accret_med'] 
    z = np.round(acc['redshifts'], 3)

    if width is not None:
        if width =='std':
            acc_width = acc['accret_std']
        elif width == 'percentile':
            lowp_array  = acc['accret_low']
            highp_array = acc['accret_high']
            acc_width = highp_array - lowp_array
        return mass_cuts, z, accret_med, acc_width
    else:
        return mass_cuts, z, accret_med 
    
def plot_data(dir, snaps, acc, axs, cmap, norm, feat=None):
    import numpy as np
    import seaborn as sns
    from scipy.stats import pearsonr
    
    """
    Params:
        dir:    str
        snaps:  list
        acc:    list 
                [acc_z:         (M,)
                 acc_mass_cuts: (N,) 
                 acc_data:      (M, N,)
                 acc_data_err:  (M, N,)]    
        feat:   str ('depth' or 'width')
        """
    tot_x, tot_xerr = [], []
    tot_z = []
    tot_y, tot_y_min, tot_y_max = [], [], []
    
    acc_z, acc_mass_cuts, acc_data, acc_err = acc[0], acc[1], acc[2], acc[3]

    for snap in snaps:

        # Load data
        data = np.load(dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
        z = np.round(data['z'], 3)
        print(f'snap {snap}: z = ', z)
        mass_cuts = data['mass_bins']
        print(mass_cuts)
        data = data['final_results']
        
        # In accretion rate, find the index corresponding to the current snap
        acc_snap_idx = np.where(acc_z == z)[0][0]
        print('acc z =', acc_z[acc_snap_idx], 'at idx', acc_snap_idx, acc_data.shape)
        
        # In accretion rate, find the index corresponding to the current mass cut
        comm_mass_cuts = np.intersect1d(acc_mass_cuts, mass_cuts)
        print(comm_mass_cuts)
        x_mass_idx = [np.where(acc_mass_cuts == m)[0][0] for m in comm_mass_cuts]
        print(x_mass_idx, acc_mass_cuts[x_mass_idx])
        y_mass_idx = [np.where(mass_cuts == m)[0][0] for m in comm_mass_cuts]
        print(y_mass_idx, mass_cuts[y_mass_idx])
        
        x = acc_data[acc_snap_idx, x_mass_idx]
        xerr = acc_err[acc_snap_idx, x_mass_idx]/2
        
        if feat == 'depth':
            y = data[y_mass_idx, 1, 1]
            y_min, y_max = data[y_mass_idx, 1, 0], data[y_mass_idx, 1, 2]
        elif feat == 'width':
            y = data[y_mass_idx, 2, 1]
            y_min, y_max = data[y_mass_idx, 2, 0], data[y_mass_idx, 2, 2]

        mask = x!=0
        x, xerr, y, y_min, y_max = x[mask], xerr[mask], y[mask], y_min[mask], y_max[mask]
        
        # Sort according to x
        x = x[np.argsort(x)]
        xerr = xerr[np.argsort(x)]
        y = y[np.argsort(x)]
        y_min = y_min[np.argsort(x)]
        y_max = y_max[np.argsort(x)]
        
        if len(x) != 0:
            
            print('x:', x)
            print('y:', y)
            print('')
            # axs.plot(x, y, color=cmap(norm(np.round(z, 2))))
            axs.errorbar(x, y, yerr=[y-y_min, y_max-y], # xerr=xerr, 
                         color=cmap(norm(np.round(z, 2))), fmt='.')
        
            tot_x = np.concatenate((tot_x, x))
            tot_xerr = np.concatenate((tot_xerr, xerr))
            tot_z.append([z]*len(x))
            
            tot_y = np.concatenate((tot_y, y))
            tot_y_min = np.concatenate((tot_y_min, y-y_min))
            tot_y_max = np.concatenate((tot_y_max, y_max-y))
            
    tot_z = np.concatenate(tot_z)
    return tot_x, tot_xerr, tot_z, tot_y, [tot_y_min, tot_y_max]

def delta_c(z):
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

# def poly_fit(X, Y, Yerr, redshifts, axs, cmap, norm, DoF):
#     """Params:
#         X: 1d array
#         Y: 1d array
#         redshifts: 1d array"""

#     for z in np.unique(redshifts):
        
#         fit_X = X[np.where(np.round(redshifts, 1) == np.round(z, 1))]
#         fit_Y = Y[np.where(np.round(redshifts, 1) == np.round(z, 1))]
#         print(fit_X)
        
#         if len(fit_X) > DoF+1:
#             # fit            
#             popt, _ = np.polyfit(fit_X, fit_Y, deg=DoF, cov=True)
#             poly = np.poly1d(popt)
#             Y_fit = poly(X)
            
#             # find the reduced chi-square
#             red_chi2 = np.sum((Y - Y_fit)**2 / Yerr**2) / (len(Y) - len(popt))
#             print(f'z = {z}, reduced chi2 = {red_chi2}')
            
#             # Plot
#             axs.plot(np.linspace(fit_X[0], fit_X[-1], 100), 
#                      poly(np.linspace(fit_X[0], fit_X[-1], 100)), 
#                      color=cmap(norm(np.round(z, 3))),
#                      label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}')
            
#             # Print the fitted polynomial
#             print(f'z = {z}:')
#             print(poly)
#     return axs