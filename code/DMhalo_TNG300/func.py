def compt_density_profile(coordinates, haloCM, halo_R_Mean200, dimension,
                          radius_range=[0.01, 5], number_of_bins=85):
    """This function computes the density profile of halos.
    
    INPUT:
    coordinates:    2D array with shape (N, 3)              [free unit]
    haloCM:         1D array with shape (3,)                [free unit]
    halo_R_Mean200: float
    dimension:      float                                   [ckpc^3/h]
    radius_range:   list with two fractional radius         [dimensionless]
    number_of_bins: float                       
    
    RETURN:
    density_bins:   1D array with shape (number of bins,)   [input unit^(-3)]
    radius_bins:    1D array with shape (number of bins,)   [input unit]
    """
    
    import numpy as np
    # radial_coordinates = np.linalg.norm(coordinates - haloCM, axis=1) 
    radial_coordinates = distance(haloCM, coordinates, dimension)
    
    # Create bins
    radius_bins = np.logspace(np.log10(radius_range[0]*halo_R_Mean200), 
                                 np.log10(radius_range[1]*halo_R_Mean200), number_of_bins)
    density_bins = np.empty((number_of_bins))
    
    # Compute the density profile 
    for i, r in enumerate(radius_bins):
        # mass enclosed in bins
        if i == 0:
            mass_bin = len(radial_coordinates[(radial_coordinates <= r)]) 
            volume_bin =  4/3 * np.pi * r**3 
        else:
            mass_bin = len(radial_coordinates[(radius_bins[i-1] < radial_coordinates) & 
                                                    (radial_coordinates <= r)]) 
            volume_bin =  4/3 * np.pi * (r**3 - radius_bins[i-1]**3) 
        # Compute density
        density_bins[i] = mass_bin / volume_bin
    
    return density_bins, radius_bins

def distance(x0, x1, dimensions):
    '''
    find distance between points in a periodic box
    x0,x1=coordinates or set of coordinates
    delta = distance between two points assuming no periodicity
    dimensions = box size *7500 ckpc/h*
    returns delta if delta<box size, else returns 1-delta as the distance between points
    '''
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def compt_density_profile_hist(coordinates, mass_weights, haloPos, radial_bins):
    
    radii = np.sqrt(np.sum(coordinates-haloPos, axis=1)**2) # [ckpc/h]
    radial_volumes = 4/3*np.pi * (radial_bins[1:]**3 - radial_bins[:-1]**3) # [(ckpc/h)^3]
    densities = np.histogram(radii,radial_bins,weights=mass_weights)[0] / radial_volumes # [Msun/h / (ckpc/h)^3]
    return densities

def select_halos(file_list, bin_start, bin_end):
    
    import numpy as np
    
    sub_list = []
    for file in file_list:
        # Load data
        data = np.load(file, allow_pickle=True).item()
        halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 MSun/h]
        
        # Apply the mass criteria 
        if (halo_M_Mean200 >= 10**bin_start) & (halo_M_Mean200 < 10**(bin_end+0.5)):
            sub_list.append(file)
    return sub_list
    
def count_halos(root_dir, fnames, bin_start, bin_end):
    
    import os
    import numpy as np
    # import illustris_python as il
    
    counts = []
    
    bin_width = 0.5
    num_bins = int((bin_end-bin_start)/bin_width)
    bins = np.arange(bin_start, bin_end, bin_width)

    for i in range(num_bins):
        count = 0
        for fname in fnames:
            # Load data
            data = np.load(os.path.join(root_dir, fname), allow_pickle=True).item()
            halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 MSun/h]
            
            # Apply the mass criteria 
            if (halo_M_Mean200 >= 10**bins[i]) & (halo_M_Mean200 < 10**(bins[i]+0.5)):
                count += 1
        counts.append(count)
    return counts
    
def stacked_density_profile(root_dir, fnames, mass_criteria, use_bootstrap=True):
    """The function computes the stacked density profiles.
    
    INPUT:
    
    root_dir:       str
    fnames:         list of filepaths without subfolders!
    mass_criteria:  list [a, b]                 [10^(10+a) MSun/h]
        
    RETURN:
    (   
    radial_centers: 1D array with shape (N,)    [dimensionless]
    medians:        1D array with shape (N,)    [dimensionless]
    errors:         2D array with shape (N, 2)  [dimensionless]
    num_halo:       int
    R200_median:    float                       [ckpc/h]
    )
    """
    
    import os
    import sys
    from tqdm import tqdm
    import numpy as np
    import matplotlib.pyplot as plt
    from unyt import G, second, megaparsec, km
    
    # Initialize the output
    density_profiles = [] # [dimensionless]
    radii = []            # [dimensionless]
    R200 = []             # [ckpc/h]
    
    ### Load all density profile data ###
    # Iterate over halos
    for fname in tqdm(fnames):
        # Load data
        data = np.load(os.path.join(root_dir, fname), allow_pickle=True).item()
        h = data['h'] # unit [100 * km / megaparsec / second]
        scale_factor = data['scale_factor']
        halo_M_Mean200 = data['halo_M_Mean200'] # [10^10 MSun/h]
        
        # Apply the mass criteria 
        if (halo_M_Mean200 >= 10**mass_criteria[0]) & (halo_M_Mean200 < 10**mass_criteria[1]):
            # radius profile
            radii = [item / data['halo_R_Mean200'] for item in data['radial_bins']] 
            R200.append(data['halo_R_Mean200'])
            
            # Compute the critical density
            Hubble = data['h'] * 100 * km / megaparsec / second
            rho_c = 3 * Hubble**2 / (8*np.pi*G) 
            rho_c.convert_to_units('Msun/kiloparsec**3') 
            rho_c = rho_c * scale_factor**3 / h**2 # [(MSun/h)/(ckpc/h)**3]
            
            # density profile
            density_profiles.append(data['densities'] / rho_c) 
           
        else:
            pass 
    
    # Compute median R200
    R200_median = np.median(R200) # [ckpc/h]
    del R200
    
    # Get the number of DM halos
    num_halo = np.array(density_profiles).shape[0]
    if num_halo == 0:
        print('no halos')
        # sys.exit()
    elif num_halo < 10:
        print(f'num of halos: {num_halo}, too few')
        # sys.exit()
    else:
        print(f'num of halos: {num_halo}')
        
    density_profiles = np.array(density_profiles).T # shape: (N radii, num_of_halos)
    
    ### Find the median density profile ###
    if use_bootstrap == True:
        medians, errors = [], []
        for idx, rho_data_at_r in enumerate(density_profiles):
            resampled_median_rho_data_at_r = bootstrap(rho_data_at_r, np.median)
            # plt.figure()
            # plt.hist(resampled_median_rho_data_at_r)
            # plt.savefig(f'result/hist/r_{idx}')
            # plt.close()
            median_with_error = np.percentile(resampled_median_rho_data_at_r, [16, 50, 84])
            # Append new data to the lists
            medians.append(median_with_error[1])
            errors.append([abs(median_with_error[1]-median_with_error[0]), 
                           abs(median_with_error[2]-median_with_error[1])])
        # Convert the list to array
        medians = np.array(medians) # shape: (N radii, 1)
        errors = np.array(errors) # shape: (N radii, 2)
    else:
        medians = np.percentile(density_profiles, 50, axis=1) # shape: (N radii, 1)
        errors = np.percentile(density_profiles, [16,84], axis=1).T # shape: (N radii, 2)
    
    # Compute radii centers
    radial_centers = np.array(radii)

    # check the return radial_centers and medians shape
    if radial_centers.shape != medians.shape:
        radial_centers = radial_centers[1:]
        
    return (radial_centers, medians, errors, num_halo, R200_median)

def bootstrap(x, statfunc, Nboots=32):
    import numpy as np
    
    x = np.array(x)
    
    resampled_stat = []
    for k in range(Nboots):
        index = np.random.randint(0, len(x), len(x))
        sample = x[index]
        stat_value = statfunc(sample)
        resampled_stat.append(stat_value)
    
    return resampled_stat

# Compute d log rho / d log r

def filter_profile(x, y, yerr):
    """
    x: shape (N,)
    y: shape (N,)
    yerr: shape (2, N,)"""

    for i in range(len(x)):
        if all(item != 0 for item in y[i:]):
            start_idx = i
            break

    # start_idx = next((i for i, x in enumerate(y) if x != 0), 0)
    return x[start_idx:], y[start_idx:], yerr[start_idx:]
    
def num_deriv(lgR, lgP):
    result = np.zeros(len(lgR))
    result[2:-2] = ( 1./12.*lgP[0:-4] - 2./3.*lgP[1:-3] + 2./3.*lgP[3:-1] - 1./12.*lgP[4:] ) / (lgR[3:-1] - lgR[2:-2])
    result[0]    = (lgP[1] - lgP[0]) / (lgR[1] - lgR[0])
    result[1]    = (lgP[2] - lgP[0]) / (lgR[2] - lgR[0])
    result[-1]   = (lgP[-1] - lgP[-2]) / (lgR[-1] - lgR[-2])
    result[-2]   = (lgP[-3] - lgP[-1]) / (lgR[-3] - lgR[-1])
    return result
 
def num_deriv_err_single(x,y,yerr):
    # Compute log y error
    lgy_err = abs(yerr / y)
    # Compute d log y / d log x error
    lgx, lgy = np.log10(x), np.log10(y)
    slope = num_deriv(lgx, lgy)
    err = abs(slope * lgy_err / lgy)
    return err

def num_deriv_err(x,y,yerrs):
    errs = np.empty((x.shape[0], 2))
    errs[:,0] = num_deriv_err_single(x,y, yerrs[:,0])
    errs[:,1] = num_deriv_err_single(x,y, yerrs[:,1])
    return errs

def float_to_str(bin_start, bin_end):
    if str(bin_start).endswith('0'):
        start = int(bin_start)
        enda, endb = str(bin_end).split('.')
        end = f'{enda}-{endb}'
    else:
        starta, startb = str(bin_start).split('.')
        start = f'{starta}-{startb}'
        end = int(bin_end)
    return start, end

def float_to_int(bin_start, bin_end):
    start = int(bin_start * 10)
    end = int(bin_end * 10)
    return start, end

def plot_profile(R200_median, radius, rho, rho_err, slope, slope_err, 
                 fitted_radius, fitted_rho, fitted_slope, 
                 mass_cut, num_halo, snap, save_dir, fname,
                 save_data = False):
    
    import os
    from matplotlib import pyplot as plt  
    plt.style.use('code/style.mplstyle')
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 15))
    axs[0].scatter(radius, rho, s=1, # color='b', 
                   label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot/h$".format(mass_cut[0]+10, mass_cut[1]+10)
                   # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                   )
    axs[0].fill_between(radius, rho-rho_err[:,0], rho+rho_err[:,1], alpha = 0.2, # color = 'b',
                        # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                        )
    axs[0].plot(fitted_radius, fitted_rho, lw=0.5, # color='salmon',
                # label=f'Fit: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                )
    
    # Plot the fitted gradients
    axs[1].scatter(radius, slope, s=1, # color='b',
                   label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot/h$".format(mass_cut[0]+10, mass_cut[1]+10)
                   # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                   )
    axs[1].fill_between(radius, slope-slope_err[:,0], slope+slope_err[:,1], alpha = 0.2, # color = 'b',
                        # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                        )
    axs[1].plot(fitted_radius, fitted_slope, lw=0.5, # color='salmon',
                # label=f'Theory: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
                )
    

    # General settings
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rho$/$\rho_c$")
    axs[0].legend()
    axs[0].set_title(f'Stacked density profiles')

    axs[1].set_xscale('log')
    axs[1].set_xlabel(r"r/$R_{200}$")
    axs[1].set_ylabel("Slope")
    axs[1].set_ylim(-6,-0)
    axs[1].legend()
    axs[1].set_title(f'Finding splashback radius at snap {snap}')
    
    # Save the plot
    start, _ = float_to_str(mass_cut[0], mass_cut[1])
    plt_dir = save_dir + f'/plot/mass_cut_{start}'
    if not os.path.exists(plt_dir):
        os.makedirs(plt_dir)
    plt.savefig(os.path.join(plt_dir, fname))
    plt.close()
    
    # Save data
    if save_data == True:
        data = {'radius': radius, 'rho': rho, 'rho_err': rho_err,
                'slope': slope, 'slope_err': slope_err, 
                'fitted_radius': fitted_radius, 'fitted_rho': fitted_rho, 'fitted_slope': fitted_slope,
                'R200_median': R200_median}
        
        data_dir = save_dir + f'/data/mass_cut_{start}'
        if not os.path.exists(data_dir):
           os.makedirs(data_dir)
        np.save(os.path.join(data_dir, fname), data)



### External functions ###

import numpy as np
from scipy.optimize import curve_fit

evaluate_profile_at_edges = np.logspace(-2, np.log10(5), 1024)
evaluate_profile_at = 0.5 * (evaluate_profile_at_edges[1:] + evaluate_profile_at_edges[:-1])

def NSW_profile(r,
                rho_0, 
                R_s):
    return rho_0 / ((1 + (r/R_s)**2)*(r/R_s))
    
def density_profile_inner(
        r,
        rho_s,
        r_s,
        alpha,
    ):  
        
        return np.log10(rho_s * np.exp(-(2.0 / alpha) * (np.power(r / r_s, alpha) - 1.0)))

def density_profile_outer(
        r,
        rho_g,
        b_e,
        S_e,
        R_200_mean,
    ):
        return np.log10(rho_g * (b_e * np.power(r / (5.0 * R_200_mean), -S_e) + 1.0))

def density_profile(
    r: float,
    rho_s: float,
    r_s: float,
    alpha: float,
    r_t: float,
    beta: float,
    gamma: float,
    rho_g: float,
    b_e: float,
    S_e: float,
    R_200_mean: float,
) -> float:

    rho_inner = rho_s * np.exp(-(2.0 / alpha) * (np.power(r / r_s, alpha) - 1.0))
    f_trans = np.power(1.0 + np.power(r / r_t, beta), -gamma / beta)
    rho_outer = rho_g * (b_e * np.power(r / (5.0 * R_200_mean), -S_e) + 1.0)

    rho = rho_inner * f_trans + rho_outer

    return np.log10(rho)

def density_gradient_profile(
    r: float,
    rho_s: float,
    r_s: float,
    alpha: float,
    r_t: float,
    beta: float,
    gamma: float,
    rho_g: float,
    b_e: float,
    S_e: float,
    R_200_mean: float,
) -> float:

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

def fit_NSW_profile(bin_centers, densities, R_200_mean):
    
    def wrapped_NSW_profile(r:float, 
                        rho_0:float, 
                        R_s:float
                        ):
        return NSW_profile(r=r, rho_0=rho_0, R_s=R_s)
    
    r = bin_centers * R_200_mean
    rho = densities / (R_200_mean ** 3)
    
    base_p0 = (max(rho), R_200_mean)
    print(base_p0)
    popt, pcov = curve_fit(NSW_profile,
                            r,
                            rho,
                            p0=base_p0,
                            maxfev=100000)
    perr = np.diag(pcov) ** 0.5
    
    def chi_square(p):
        return np.sum((p - densities)**2) / (len(p) - len(popt))

    # Did we actually get a good fit? If not, we should dump this bootstrapping.
    predicted_values = wrapped_NSW_profile(bin_centers * R_200_mean, *popt)
    new_chi_square = chi_square(predicted_values)
    # If we get a worse fit after tuning, cancel this one
    old_chi_square = chi_square(wrapped_NSW_profile(bin_centers * R_200_mean, *base_p0))

    # if old_chi_square < new_chi_square:
    #     raise RuntimeError("Extremely poor fit for this bootstrap")
    # elif new_chi_square > 2.0:
    #     print("Bad Fit!")
    #     raise RuntimeError("Extremely poor fit for this bootstrap")

    return (
        popt,perr,
        )
    
def fit_gradient_parametric(bin_centers, densities, density_errors, gradient, gradient_errors, R_200_mean):
    """
    Fits a profile in log-log spce based on Equation 6 in
    O'Neil et al. 2021.

    First fits the inner and outer profile separately, then freezes those
    parameters, fitting the transfer function. Finally, all parameters are
    released to fine-tune the fit.

    Doesn't make use of sigma.
    """

    log_rho = np.log10(densities / (R_200_mean ** 3))

    log_rho_upper = np.log10((densities + density_errors) / (R_200_mean ** 3))
    log_rho_lower = np.log10(np.maximum((densities - density_errors) / (R_200_mean ** 3), 0.01 * densities / (R_200_mean ** 3)))

    log_rho_error = 0.5 * (log_rho_upper - log_rho_lower)

    def wrapped_profile_gradient(
        r: float,
        rho_s: float,
        r_s: float,
        alpha: float,
        r_t: float,
        beta: float,
        gamma: float,
        rho_g: float,
        b_e: float,
        S_e: float,
    ):
        return density_gradient_profile(
            r=r,
            rho_s=rho_s,
            r_s=r_s,
            r_t=r_t,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho_g=rho_g,
            b_e=b_e,
            S_e=S_e,
            R_200_mean=R_200_mean,
        )

    def wrapped_profile(
        r: float,
        rho_s: float,
        r_s: float,
        alpha: float,
        r_t: float,
        beta: float,
        gamma: float,
        rho_g: float,
        b_e: float,
        S_e: float,
    ):
        return density_profile(
            r=r,
            rho_s=rho_s,
            r_s=r_s,
            r_t=r_t,
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            rho_g=rho_g,
            b_e=b_e,
            S_e=S_e,
            R_200_mean=R_200_mean,
        )

    # First, fit inner profile.
    r = bin_centers * R_200_mean

    # Global mask
    global_mask = bin_centers > 0.02

    # First, fit inner profile.
    r = bin_centers * R_200_mean

    inner_mask = np.logical_and(bin_centers > 0.02, bin_centers < 0.8)

    popt_inner, pcov_inner = curve_fit(
        density_profile_inner,
        r[inner_mask],
        log_rho[inner_mask],
        p0=(10**log_rho[0], R_200_mean, 0.5),
        sigma=log_rho_error[inner_mask],
        # For some reason bounds make this go very wrong.
        # bounds=(
        #     [1e-10 * 10**log_rho[0], 0.001 * R_200_mean, 0.0],
        #     [np.inf, 10.0 * R_200_mean, 1]
        # ),
        maxfev=100000,
    )

    def wrapped_outer(
        r: float,
        rho_g,
        b_e,
        S_e,
    ):
        return density_profile_outer(
                r=r,
                rho_g=rho_g,
                b_e=b_e,
                S_e=S_e,
                R_200_mean=R_200_mean,
            )

    outer_mask = bin_centers > 2.0

    # Using bounds here messes this up because it can no longer
    # use lm, and instead uses trf, unless they are very tight.
    # In particular, our requirement that b_e > 1.0 is required.
    popt_outer, pcov_outer = curve_fit(
        wrapped_outer,
        r[outer_mask],
        log_rho[outer_mask],
        p0=(10 ** log_rho[-1],
            2.0,
            2.0,),
        sigma=log_rho_error[outer_mask],
        bounds=(
            [0.01 * 10 ** log_rho[-1], 0.1, 0.0],
            [10 * 10 ** log_rho[-1], 5.0, 5.0]
        ),
        maxfev=100000,
    )

    # return (
    #     evaluate_profile_at,(
    #     10 ** wrapped_outer(evaluate_profile_at * R_200_mean, *popt_outer)
    #     #+ 10 ** density_profile_inner(evaluate_profile_at * R_200_mean, *popt_inner)
    #     )
    #     * R_200_mean ** 3,
    # )

    # Now fit f_trans separately

    def wrapped_ftrans_only(
        r: float,
        r_t: float,
        beta: float,
        gamma: float,
    ):
        return wrapped_profile(
            r,
            *popt_inner,
            r_t,
            beta,
            gamma,
            *popt_outer,
        )

    middle_mask = np.logical_and(
        bin_centers > 0.5, bin_centers < 2.0
    )

    popt_ft, pcov_ft = curve_fit(
        wrapped_ftrans_only,
        r[middle_mask],
        log_rho[middle_mask],
        p0=(R_200_mean,
            2.0,
            4.0,),
        maxfev=100000,
        sigma=log_rho_error[middle_mask],
        bounds=(
            [0.01 * R_200_mean, 0.0, 0.0],
            [10.0 * R_200_mean, 10.0, 5.0]
        )
    )

    base_p0 = (
        *popt_inner, *popt_ft, *popt_outer,
    )

    change_frac = 4.0

    base_lower = [x / change_frac for x in base_p0]
    base_upper = [x * change_frac for x in base_p0]



    # p0 = p0_full if p0_full is not None else base_p0

    # Now fit the gradient separately.

    popt, pcov = curve_fit(
        wrapped_profile_gradient,
        bin_centers * R_200_mean,
        gradient,
        p0=base_p0,
        maxfev=100000,
        bounds=[base_lower, base_upper],
        sigma=gradient_errors,
    )

    # popt = base_p0

    def chi_square(p):
        return np.sum(((p - gradient)/ gradient_errors)**2) / (len(p) - len(popt))

    # Did we actually get a good fit? If not, we should dump this bootstrapping.
    predicted_values = wrapped_profile_gradient(bin_centers * R_200_mean, *popt)
    new_chi_square = chi_square(predicted_values)

    # If we get a worse fit after tuning, cancel this one
    old_chi_square = chi_square(wrapped_profile_gradient(bin_centers * R_200_mean, *base_p0))

    if old_chi_square < new_chi_square:
        raise RuntimeError("Extremely poor fit for this bootstrap")
    # elif new_chi_square > 2.0:
    #     print("Bad Fit!")
    #     raise RuntimeError("Extremely poor fit for this bootstrap")

    return (
        evaluate_profile_at,
        wrapped_profile_gradient(evaluate_profile_at * R_200_mean, *popt),
    )