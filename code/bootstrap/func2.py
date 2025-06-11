### External functions ###

import numpy as np
from scipy.optimize import curve_fit
from lmfit import Model, Parameters, minimize
import numpy as np

evaluate_profile_at_edges = np.logspace(-2, np.log10(5), 1024)
evaluate_profile_at = 0.5 * (evaluate_profile_at_edges[1:] + evaluate_profile_at_edges[:-1])


def density_profile_fk(
        r,
        rho_s,
        r_s,
        alpha,
        r_t,
        beta
    ):  
        
        return np.log10(rho_s * np.exp(-(2.0/alpha)*(np.power(r/r_s, alpha) - (1.0/beta)*(r/r_t)**beta)))


def gradient_profile_fk(
    r: float,
    rho_s: float,
    r_s: float,
    alpha: float,
    r_t: float,
    beta: float,) -> float:
    
    rho = rho_s * np.exp(-(2.0/alpha)*(np.power(r/r_s, alpha) - (1.0/beta)*(r/r_t)**beta))
    
    # Different functional form for inner profile
    # -------------------------------------------
    d_rho_dr = ( - (2.0/r_s)*np.power(r/r_s, alpha-1) - (1.0/r_t)*np.power(r/r_t, beta-1) ) * rho

    return r/rho * d_rho_dr
        
        
def fit_profile_parametric_fk(bin_centers, densities, R_200_mean):

    log_rho = np.log10(densities / (R_200_mean ** 3))
    r = bin_centers * R_200_mean
    
    try:
        popt1, _ = curve_fit(density_profile_fk, r, log_rho,
                                    p0=(10**log_rho[0], 
                                        R_200_mean, 0.18,
                                        R_200_mean, 0.3,), # Different functional form
                                    # sigma=log_rho_error[inner_mask],
                                    bounds=(
                                        [0.01*10**log_rho[0], 0.01*R_200_mean, 0.03, 0.5*R_200_mean, 0.1],
                                        [20*10**log_rho[0],   4.5*R_200_mean, 0.4,  3*R_200_mean, 10]
                                           ),
                                    maxfev=10000000,
                                )       
        # Fit again
        popt2, _ = curve_fit(density_profile_fk, r, log_rho, p0=popt1, maxfev=10000000)   
        
        # Compare the chi square
        chi2_1 = np.sum((density_profile_fk(r, *popt1) - log_rho) ** 2)
        chi2_2 = np.sum((density_profile_fk(r, *popt2) - log_rho) ** 2)
        if chi2_2 < chi2_1: 
            print('Successfully found optimal params! ')
            return (
                evaluate_profile_at,
                10 ** density_profile_fk(evaluate_profile_at * R_200_mean, *popt2) * R_200_mean ** 3,
                popt2,
            )  
        else:
            print('Warning: Optimal parameters not found, using the first fit.')
            return (
                evaluate_profile_at,
                10 ** density_profile_fk(evaluate_profile_at * R_200_mean, *popt1) * R_200_mean ** 3,
                popt1,
            )
                        
    # Inner curve fit error
    except RuntimeError as e:
        print(f"Warning: Optimal parameters not found. Error: {e}")
        return (
            evaluate_profile_at,
            np.array([0]),
            np.array([0, 0, 0, 0, 0])
        )       
        
        
def fit_gradient_parametric_fk(bin_centers, gradients, R_200_mean, init_p0):

    r = bin_centers * R_200_mean
    print('print p0', init_p0)
    try:
        popt, _ = curve_fit(gradient_profile_fk,
                            r,gradients,
                            p0=init_p0,
                            # sigma=log_rho_error,
                            maxfev=1000000,
                            # bounds=(
                            #     [0.01*10**init_p0[0], 0.1*init_p0[1], 0.1,  0.01*R_200_mean, 0.01],
                            #     [20*10**init_p0[0],   50*init_p0[1],  20.0, 100*R_200_mean,  50.0]
                            #         ),
                        )    
           
        print('Successfully found optimal params! ')
        return (
            evaluate_profile_at,
            gradient_profile_fk(evaluate_profile_at*R_200_mean, *popt) * R_200_mean ** 3,
            popt,
        )
        
    except RuntimeError as e:
        print(f"Warning: Optimal parameters not found. Error: {e}")
        return (
            evaluate_profile_at,
            np.full((1023), 0),
            np.array([0, 0, 0, 0, 0,]),
        ) 
        
        
        
# def density_profile_fk(params, r, log_rho): 
#     rho_s = params['rho_s'].value
#     r_s = params['r_s'].value
#     alpha = params['alpha'].value
#     r_t = params['r_t'].value
#     beta = params['beta'].value
#     return np.abs(log_rho - (np.log10(rho_s * np.exp(-(2.0/alpha)*(np.power(r/r_s, alpha) - (1.0/beta)*(r/r_t)**beta)))))
    
# def fit_profile_parametric_fk(bin_centers, densities, R_200_mean):
#     log_rho = np.log10(densities / (R_200_mean ** 3))
#     r = bin_centers * R_200_mean

#     params = Parameters()
#     init_rho_s = 10 ** log_rho[0]
#     params.add('rho_s', value=init_rho_s, min=0.01*init_rho_s, max=20*init_rho_s)
#     params.add('r_s',   value=R_200_mean, min=0.01*R_200_mean, max=4.5*R_200_mean)
#     params.add('alpha', value=0.18,       min=0.03,            max=0.4)
#     params.add('r_t',   value=R_200_mean, min=0.5*R_200_mean,  max=3*R_200_mean)
#     params.add('beta',  value=0.3,        min=0.1,             max=10)
    
#     result = minimize(density_profile_fk, params, args=(r, log_rho))
#     print(result)
#     if result.success:
#         best_vals = result.best_values
#         y_fit = 10 ** density_profile_fk(evaluate_profile_at * R_200_mean, **best_vals) * R_200_mean ** 3
#         return evaluate_profile_at, y_fit, best_vals
#     else:
#         return evaluate_profile_at, np.array([0]), {k: 0 for k in params.keys()}
    
# def fit_gradient_parametric_fk(bin_centers, gradients, R_200_mean, init_p0):
    r = bin_centers * R_200_mean
    model = Model(gradient_profile_fk)

    params = Parameters()
    for name, val in zip(['rho_s', 'r_s', 'alpha', 'r_t', 'beta'], init_p0):
        params.add(name, value=val)

    result = model.fit(gradients, params, r=r)
    if result.success:
        best_vals = result.best_values
        y_fit = gradient_profile_fk(evaluate_profile_at * R_200_mean, **best_vals) * R_200_mean ** 3
        return evaluate_profile_at, y_fit, best_vals
    else:
        return evaluate_profile_at, np.zeros_like(evaluate_profile_at), {k: 0 for k in params.keys()}