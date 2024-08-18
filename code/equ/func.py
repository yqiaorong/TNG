def log_error(y, yerr):
    import numpy as np
    
    log_y_upper = np.log10((y + yerr))
    log_y_lower = np.log10(np.maximum(y - yerr), 0.01*y)

    return 0.5 * (log_y_upper - log_y_lower)

def partial_func(x, a, b, c):
    return a*x**2 + b*x + c

def model_func(x, a, b, c, d, e, f):
    logM, z = x
    return a*logM**2 + b*logM + c + d*z**2 + e*z + f

def fit_profile_parametric(mass, z, y, yerr, p0):
    """Fit profile in log space of mass. """  

    def wrapped_profile(a: float,
                        b: float,
                        c: float,
                        d: float,
                        e: float,
                        f: float,):
        return density_profile(a=a, b=b, c=c, d=d, e=e, f=f)
    
    # Stack two inputs
    x = np.vstack((np.log10(mass), z))
    
    popt, _ = curve_fit(model_func, x, y, p0=p0,
                        sigma=log_rho_error[inner_mask],
                        # maxfev=100000,
                        )


    change_frac = 1.1

    # base_lower = [x / change_frac if x >=0 else x * change_frac for x in base_p0]
    # base_upper = [x * change_frac if x >=0 else x / change_frac for x in base_p0]
    
    base_lower = [x / change_frac if x > 0 else (-0.1 if x == 0 else x * change_frac) for x in base_p0]
    base_upper = [x * change_frac if x > 0 else (0.1 if x == 0 else x / change_frac) for x in base_p0]
    # print(base_lower)
    # print(base_upper)

    # p0 = p0_full if p0_full is not None else base_p0


    popt, _ = curve_fit(
        wrapped_profile,
        bin_centers[global_mask] * R_200_mean,
        log_rho[global_mask],
        p0=base_p0,
        maxfev=100000,
        bounds=[base_lower, base_upper],
        # sigma=log_rho_error[global_mask],
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