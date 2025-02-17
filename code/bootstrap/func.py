def count_halos(array_1d, bin_start_list, bin_end_list):
    import numpy as np
    
    counts_list = []
    for bin_start, bin_end in zip(bin_start_list, bin_end_list):
        count = int(np.sum((array_1d >= bin_start) 
                         & (array_1d < bin_end)))
        counts_list.append(count)
    return counts_list

def stacked_density_profile(radii_2d_array, densities_2d_array, target_1d_array, R200_1d_array, 
                            bin_start, bin_end, rho_c, use_bootstrap=True):
    """The function computes the stacked density profiles.
    
    INPUT:
    
    radii_2d_array:     array (N_halos, N_radial_bins)       [kpc]
    densities_2d_array: array (N_halos, N_radial_bins)       [Msun / (kpc)^3]
    M200_1d_array:      array (N_halos,)                     [10^10 MSun]
    R200_1d_array:      array (N_halos,)                     [kpc]
    mass_bin_start:     float                                [10^(10+a) MSun]
    h:                  float                                  
    scale_factor:       float
    rho_c:              float                                [Msun / (kpc)^3]
    
    RETURN:
    (   
    radial_centers:     array (N_radial_bins)                [dimless // r200]
    medians:            array (N_radial_bins)                [dimless // rho_c]
    errors:             array with shape (N_radial_bins, 2)  [dimless // rho_c]
    num_halo:           int
    R200_median:        float                                [kpc]
    )
    """
    
    import numpy as np
    # from unyt import G, second, megaparsec, km
    
    # Initialize the output
    densities_list = [] # [dimensionless]
    radii = []          # [dimensionless]
    radius200_list = [] # [kpc]
    
    ### Load all density profile data ###
    # Iterate over halos
    for radii, densities, target, radius200 in zip(radii_2d_array, densities_2d_array, target_1d_array, R200_1d_array):
        # Load data

        # Apply the mass criteria 
        if (target >= bin_start) & (target < bin_end):
            # radius profile
            radii_dimless = radii / radius200 # [dimensionless]
            radius200_list.append(radius200)
            
            # # Compute the critical density
            # Hubble = h * 100 * km / megaparsec / second
            # rho_c = 3 * Hubble**2 / (8*np.pi*G) 
            # rho_c.convert_to_units('Msun/kiloparsec**3') 
            # rho_c = rho_c * scale_factor**3 / h**2   # [(MSun/h)/(ckpc/h)**3]
            
            # density profile
            densities_list.append(densities / rho_c)   # [dimensionless]
        else:
            pass 
    
    # Compute median R200
    radius200_median = np.median(radius200_list) # [kpc]
    del radius200_list
    
    # Get the number of DM halos
    num_halo = np.array(densities_list).shape[0]
    if num_halo == 0:
        print('no halos')
        # sys.exit()
    elif num_halo < 10:
        print(f'num of halos: {num_halo}, too few')
        # sys.exit()
    else:
        print(f'num of halos: {num_halo}')
        
    select_densities_array = np.array(densities_list).T # (N_radial_bins, N_halos)
    del densities_list
    
    ### Find the median density profile ###
    if use_bootstrap == True:
        medians, errors = [], []
        for rho_data_at_r in select_densities_array:
            resampled_median_rho_data_at_r = bootstrap(rho_data_at_r, np.median)
            median_with_error = np.percentile(resampled_median_rho_data_at_r, [16, 50, 84])
            # Append new data to the lists
            medians.append(median_with_error[1])
            errors.append([abs(median_with_error[1]-median_with_error[0]),  # lower bound error
                           abs(median_with_error[2]-median_with_error[1])]) # upper bound error
        # Convert the list to array
        medians = np.array(medians) # shape: (N radii, 1)
        errors = np.array(errors)   # shape: (N radii, 2)
    else:
        medians = np.percentile(densities_list, 50, axis=1)       # shape: (N radii, 1)
        errors = np.percentile(densities_list, [16,84], axis=1).T # shape: (N radii, 2)
    
    # Compute radii centers
    radii_dimless = np.array(radii_dimless)

    # check the return radial_centers and medians shape
    if radii_dimless.shape != medians.shape:
        radii_dimless = radii_dimless[1:]
        
    return (radii_dimless, medians, errors, num_halo, radius200_median)

def compute_median(target_1d_array, bin_start, bin_end):
    import numpy as np
   
    # Select masses in the bin
    indices = np.where((target_1d_array >= bin_start) & (target_1d_array < bin_end))
    values = target_1d_array[indices]
    # Compute the median mass
    median = np.median(values)
    return median
    
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
    import numpy as np
    
    result = np.zeros(len(lgR))
    result[2:-2] = ( 1./12.*lgP[0:-4] - 2./3.*lgP[1:-3] + 2./3.*lgP[3:-1] - 1./12.*lgP[4:] ) / (lgR[3:-1] - lgR[2:-2])
    result[0]    = (lgP[1] - lgP[0]) / (lgR[1] - lgR[0])
    result[1]    = (lgP[2] - lgP[0]) / (lgR[2] - lgR[0])
    result[-1]   = (lgP[-1] - lgP[-2]) / (lgR[-1] - lgR[-2])
    result[-2]   = (lgP[-3] - lgP[-1]) / (lgR[-3] - lgR[-1])
    return result
 
def num_deriv_err_single(x,y,yerr):
    import numpy as np
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

# def plot_profile(R200_median, radius, rho, rho_err, slope, slope_err, 
#                  fitted_radius, fitted_rho, fitted_slope, 
#                  mass_cut, num_halo, snap, fname,
#                  save_dir=None, save_data = False):
    
#     import os
    
#     # from matplotlib import pyplot as plt  
#     # plt.style.use('code/style.mplstyle')
    
#     # fig, axs = plt.subplots(2, 1, figsize=(10, 15))
#     # axs[0].scatter(radius, rho, s=1, # color='b', 
#     #                label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot$".format(mass_cut[0]+10, mass_cut[1]+10)
#     #                # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #                )
#     # axs[0].fill_between(radius, rho-rho_err[:,0], rho+rho_err[:,1], alpha = 0.2, # color = 'b',
#     #                     # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #                     )
#     # axs[0].plot(fitted_radius, fitted_rho, lw=0.5, # color='salmon',
#     #             # label=f'Fit: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #             )
    
#     # # Plot the fitted gradients
#     # axs[1].scatter(radius, slope, s=1, # color='b',
#     #                label=r"mass = $10^{{{:.1f}}}$ ~ $10^{{{:.1f}}}$ $M_\odot$".format(mass_cut[0]+10, mass_cut[1]+10)
#     #                # label=f'Data: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #                )
#     # axs[1].fill_between(radius, slope-slope_err[:,0], slope+slope_err[:,1], alpha = 0.2, # color = 'b',
#     #                     # label=f'Errorbar: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #                     )
#     # axs[1].plot(fitted_radius, fitted_slope, lw=0.5, # color='salmon',
#     #             # label=f'Theory: mass bin 10^{mass_cut[0]+10} ~ 10^{mass_cut[1]+10} Msun/h: {num_halo} halos'
#     #             )
    

#     # # General settings
#     # axs[0].set_xscale('log')
#     # axs[0].set_yscale('log')
#     # axs[0].set_ylabel(r"$\rho$/$\rho_c$")
#     # axs[0].legend()
#     # axs[0].set_title(f'Stacked density profiles')

#     # axs[1].set_xscale('log')
#     # axs[1].set_xlabel(r"r/$R_{200}$")
#     # axs[1].set_ylabel("Slope")
#     # axs[1].set_ylim(-6,-0)
#     # axs[1].legend()
#     # axs[1].set_title(f'Finding splashback radius at snap {snap}')
    
#     # Save the plot
#     start, _ = float_to_str(mass_cut[0], mass_cut[1])
#     # plt_dir = save_dir + f'/plot/mass_cut_{start}'
#     # if not os.path.exists(plt_dir):
#     #     os.makedirs(plt_dir)
#     # plt.savefig(os.path.join(plt_dir, fname))
#     # plt.close()
    
#     # Save data
#     if save_data == True:
#         data = {'radius': radius, 'rho': rho, 'rho_err': rho_err,
#                 'slope': slope, 'slope_err': slope_err, 
#                 'fitted_radius': fitted_radius, 'fitted_rho': fitted_rho, 'fitted_slope': fitted_slope,
#                 'R200_median': R200_median}
        
#         data_dir = save_dir + f'/data/mass_cut_{start}'
#         if not os.path.exists(data_dir):
#            os.makedirs(data_dir)
#         np.save(os.path.join(data_dir, fname), data)


### External functions ###

import numpy as np
from scipy.optimize import curve_fit

evaluate_profile_at_edges = np.logspace(-2, np.log10(5), 1024)
evaluate_profile_at = 0.5 * (evaluate_profile_at_edges[1:] + evaluate_profile_at_edges[:-1])

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

def fit_profile_parametric(bin_centers, densities, density_errors, R_200_mean):
    """
    Fits a profile in log-log spce based on Equation 6 in
    O'Neil et al. 2021.

    First fits the inner and outer profile separately, then freezes those
    parameters, fitting the transfer function. Finally, all parameters are
    released to fine-tune the fit.

    Doesn't make use of sigma.
    """
       
    global p0_full

    log_rho = np.log10(densities / (R_200_mean ** 3))

    log_rho_upper = np.log10((densities + density_errors) / (R_200_mean ** 3))
    log_rho_lower = np.log10(np.maximum((densities - density_errors) / (R_200_mean ** 3), 0.01 * densities / (R_200_mean ** 3)))

    log_rho_error = 0.5 * (log_rho_upper - log_rho_lower)

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


    # Global mask
    global_mask = bin_centers > 0.02

    ### First, fit inner profile. ###
    r = bin_centers * R_200_mean

    inner_mask = np.logical_and(bin_centers > 0.02, bin_centers < 0.8)
    
    try:
        popt_inner, _ = curve_fit(
            density_profile_inner,
            r[inner_mask],
            log_rho[inner_mask],
            p0=(10**log_rho[0], R_200_mean, 1),
            # sigma=log_rho_error[inner_mask],
            # # For some reason bounds make this go very wrong.
            # bounds=(
            #     [1e-10 * 10**log_rho[0], 0.001 * R_200_mean, 0.0],
            #     [np.inf, 10.0 * R_200_mean, 1]
            # ),
            maxfev=100000,
        )       
        
        ### Second, fit outer profile. ###
        
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
        
        try:
            popt_outer, _ = curve_fit(
                wrapped_outer,
                r[outer_mask],
                log_rho[outer_mask],
                p0=(10 ** log_rho[-1],
                    2.0,
                    2.0,),
                # sigma=log_rho_error[outer_mask],
                bounds=(
                    [0.01 * 10 ** log_rho[-1], 1.0, 1.0],
                    [10 * 10 ** log_rho[-1], 5.0, 5.0]
                    ),
                maxfev=100000,
            )

            ### Now fit f_trans separately ###
            
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

            middle_mask = np.logical_and(bin_centers > 0.5, bin_centers < 2.0)
            
            try:
                popt_ft, _ = curve_fit(
                    wrapped_ftrans_only,
                    r[middle_mask],
                    log_rho[middle_mask],
                    p0=(R_200_mean,
                        2,
                        4,),
                    maxfev=100000,
                    # sigma=log_rho_error[middle_mask],
                    bounds=(
                        [0.1 * R_200_mean, 1.0, 1.0],
                        [2.0 * R_200_mean, 5.0, 12.0]
                    )
                )

                base_p0 = (
                    *popt_inner, *popt_ft, *popt_outer,
                )

                change_frac = 1.1

                # base_lower = [x / change_frac if x >=0 else x * change_frac for x in base_p0]
                # base_upper = [x * change_frac if x >=0 else x / change_frac for x in base_p0]
                
                base_lower = [x / change_frac if x > 0 else (-0.1 if x == 0 else x * change_frac) for x in base_p0]
                base_upper = [x * change_frac if x > 0 else (0.1 if x == 0 else x / change_frac) for x in base_p0]
                # print(base_lower)
                # print(base_upper)

                # p0 = p0_full if p0_full is not None else base_p0
                
                ### Fit the entire profile ###
                try:
                    popt, _ = curve_fit(
                        wrapped_profile,
                        bin_centers[global_mask] * R_200_mean,
                        log_rho[global_mask],
                        p0=base_p0,
                        maxfev=100000,
                        bounds=[base_lower, base_upper],
                        # sigma=log_rho_error[global_mask],
                    )

                    # def chi_square(p):
                    #     return np.sum(((p[global_mask] - log_rho[global_mask])/ log_rho_error[global_mask])**2) / (len(p[global_mask]) - len(popt))

                    # # Did we actually get a good fit? If not, we should dump this bootstrapping.
                    # predicted_values = wrapped_profile(bin_centers * R_200_mean, *popt)
                    # new_chi_square = chi_square(predicted_values)

                    # # If we get a worse fit after tuning, cancel this one
                    # old_chi_square = chi_square(wrapped_profile(bin_centers * R_200_mean, *base_p0))

                    # if old_chi_square < new_chi_square:
                    #     # Just use base_p0
                    #     popt = base_p0
                    #     if old_chi_square > 2.0:
                    #         raise RuntimeError("Extremely poor fit for this bootstrap")
                    # if new_chi_square > 2.0:
                    #     raise RuntimeError("Unable to find good fit for this bootstrap")
                    
                    print('Successfully found optimal params! ')
                    return (
                        evaluate_profile_at,
                        10 ** wrapped_profile(evaluate_profile_at * R_200_mean, *popt) * R_200_mean ** 3,
                    )
                
                # Entire curve fit error
                except RuntimeError as e:
                    print(f"Warning: Optimal parameters not found for entire profile. Error: {e}")
                    return (
                        evaluate_profile_at,
                        np.array([0]),
                    )   
            
            # Ftrans curve fit error
            except RuntimeError as e:
                print(f"Warning: Optimal parameters not found for ftrans profile. Error: {e}")
                return (
                    evaluate_profile_at,
                    np.array([0]),
                )           
            
        # Outer curve fit error
        except RuntimeError as e:
            print(f"Warning: Optimal parameters not found for outer profile. Error: {e}")
            return (
                evaluate_profile_at,
                np.array([0]),
            )
        
    # Inner curve fit error
    except RuntimeError as e:
        print(f"Warning: Optimal parameters not found for inner profile. Error: {e}")
        return (
            evaluate_profile_at,
            np.array([0]),
        )