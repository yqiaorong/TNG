def compt_density_profile(coordinates, haloCM, halo_R_Mean200, 
                          radius_range=[0.01, 5], number_of_bins=85):
    import numpy as np
    radial_coordinates = np.linalg.norm(coordinates - haloCM, axis=1) 
    
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

def stacked_density_profile(file_list, mass_criteria, use_bootstrap=True):
    
    import numpy as np
    # from scipy.constants import G
    from unyt import Msun, G, second, megaparsec, km
    
    # Initialize the output
    density_profiles = []
    radii = []
    R200 = []
    ### Load all density profile data ###
    # Iterate over halos
    for file in file_list:
        # load data
        data = np.load(file, allow_pickle=True).item()
        h = data['h']
        halo_M_Mean200 = data['halo_M_Mean200'] * h
        # apply the mass criteria 
        if (halo_M_Mean200 >= 10**mass_criteria[0]) & (halo_M_Mean200 < 10**mass_criteria[1]):
            # radius profile
            radii = data['radial_bins'] / data['halo_R_Mean200']
            R200.append(data['halo_R_Mean200'])
            
            # Compute the critical density
            H = data['h'] * 100 * km / megaparsec / second
            rho_c = 3 * H**2 / (8*np.pi*G) 
            rho_c.convert_to_units('Msun/kiloparsec**3')
            
            # density profile
            density_profiles.append(data['densities']/rho_c)
           
        else:
            pass 
    
    # Compute median R200
    R200_median = np.median(R200)
    del R200
    # Get the number of DM halos
    num_halo = np.array(density_profiles).shape[0]
    if num_halo == 0:
        print('no halos')
    else:
        print(f'num of halos: {num_halo}')
        
    density_profiles = np.array(density_profiles).T # shape: (N radii, num_of_halos)
    
    ### Find the median density profile ###
    if use_bootstrap == True:
        medians, errors = [], []
        for rho_data_at_r in density_profiles:
            resampled_median_rho_data_at_r = bootstrap(rho_data_at_r, np.median)
            median_with_error = np.percentile(resampled_median_rho_data_at_r, [16, 50, 84])
            # Append new data to the lists
            medians.append(median_with_error[1])
            errors.append([median_with_error[1]-median_with_error[0], 
                           median_with_error[2]-median_with_error[1]])
        # Convert the list to array
        medians = np.array(medians) # shape: (N radii, 1)
        errors = np.array(errors) # shape: (N radii, 2)
    else:
        medians = np.percentile(density_profiles, 50, axis=1) # shape: (N radii, 1)
        errors = np.percentile(density_profiles, [16,84], axis=1).T # shape: (N radii, 2)
    
    ### Compute radii centers
    radial_centers = radii
    # radii = np.insert(radii, 0, 0)
    # radial_centers = 10**( (np.log10(radii[1:])+np.log10(radii[:-1]))/2 )
    
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

def gradient(r, rho, rho_err=None):
    import numpy as np
    
    slopes, errs = [], []
    for i in range(r.shape[0]):
        if i >= 4:
            slope = (1/12 * np.log10(rho[i-4]) - 2/3 * np.log10(rho[i-3]) + 
                2/3 * np.log10(rho[i-1]) - 1/12 * np.log10(rho[i])) / (
                    np.log10(r[i]) - np.log10(r[i-4]))
            slopes.append(slope)
            # Compute errs
            if isinstance(rho_err, np.ndarray):
                err = abs(slope) * ((rho_err[i-4]/rho[i-4])*(1/12) + 
                            (rho_err[i-3]/rho[i-3])*(2/3) + 
                            (rho_err[i-2]/rho[i-1])*(2/3) + 
                            (rho_err[i]/rho[i])*(1/12)) 
                errs.append(err)

    return r[2:-2], np.array(slopes), np.array(errs)

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
        p0=(            10 ** log_rho[-1],
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


    popt, pcov = curve_fit(
        wrapped_profile,
        bin_centers[global_mask] * R_200_mean,
        log_rho[global_mask],
        p0=base_p0,
        maxfev=100000,
        bounds=[base_lower, base_upper],
        sigma=log_rho_error[global_mask],
    )

    def chi_square(p):
        return np.sum(((p[global_mask] - log_rho[global_mask])/ log_rho_error[global_mask])**2) / (len(p[global_mask]) - len(popt))

    # Did we actually get a good fit? If not, we should dump this bootstrapping.
    predicted_values = wrapped_profile(bin_centers * R_200_mean, *popt)
    new_chi_square = chi_square(predicted_values)

    # If we get a worse fit after tuning, cancel this one
    old_chi_square = chi_square(wrapped_profile(bin_centers * R_200_mean, *base_p0))

    # if old_chi_square < new_chi_square:
    #     # Just use base_p0
    #     popt = base_p0
    #     if old_chi_square > 2.0:
    #         raise RuntimeError("Extremely poor fit for this bootstrap")
    # if new_chi_square > 2.0:
    #     raise RuntimeError("Unable to find good fit for this bootstrap")

    return (
        evaluate_profile_at,
        10 ** wrapped_profile(evaluate_profile_at * R_200_mean, *popt)
        * R_200_mean ** 3,
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
        p0=(            10 ** log_rho[-1],
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
        # bounds=[base_lower, base_upper],
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