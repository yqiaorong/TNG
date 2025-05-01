from scipy.optimize import curve_fit
import os
import numpy as np

def log_NFW_profile(r,
                rho_0, 
                R_s):
    return np.log10(rho_0 / ((1 + r/R_s)**2 * (r/R_s)))
    
def fit_log_NFW_profile(args, bin_centers, densities, R_200_mean, M_200_mean, plot=False, idx=None):
    
    def wrapped_log_NFW_profile(r:float, 
                                rho_0:float, 
                                R_s:float
                                ):
        return log_NFW_profile(r=r, rho_0=rho_0, R_s=R_s)
    
    def chi_square(p):
        return np.sum((p - densities)**2) / (len(p) - len(popt))
    
    base_p0 = (max(densities), R_200_mean)

    popt, pcov = curve_fit(log_NFW_profile,
                            bin_centers,
                            densities,
                            p0=base_p0,
                            maxfev=1000000)
    perr = np.diag(pcov) ** 0.5
    
    # if plot == True:
    
    #     # Did we actually get a good fit? If not, we should dump this bootstrapping.
    #     predicted_values = wrapped_log_NFW_profile(bin_centers, *popt)
        # new_chi_square = chi_square(predicted_values)
        # # If we get a worse fit after tuning, cancel this one
        # old_chi_square = chi_square(wrapped_log_NFW_profile(bin_centers, *base_p0))
        # print(new_chi_square, old_chi_square)
        
        # if old_chi_square < new_chi_square:
        #     pass
        #   #  raise RuntimeError("Extremely poor fit for this bootstrap")
        # elif new_chi_square > 2.0:
        #     print("Bad Fit!")
        #   #  raise RuntimeError("Extremely poor fit for this bootstrap")
        
        # # Plot densities
        # from matplotlib import pyplot as plt
        # fig, axs = plt.subplots(1,1)
        # axs.scatter(bin_centers, densities)
        # axs.plot(bin_centers, predicted_values)
        # print('bin centers', bin_centers.shape)
        # # Plot vertical line at R_s
        # axs.axvline(x=popt[1], color='red', linestyle='--', label=f'R_s = {popt[1]:.2f}')
        # axs.axvline(x=R_200_mean, color='blue', linestyle='--', label=f'R_200 = {R_200_mean:.2f}')
        # axs.text(0.1, 0.1, f'M200 = {M_200_mean:.2e}', transform=axs.transAxes)
        # axs.legend()
        # axs.set_yscale('log')
        # axs.set_xscale('log')
        # save_dir = args.sim
        # if not os.path.exists(f'output/{save_dir}'):
        #     os.makedirs(f'output/{save_dir}')
        # plt.savefig(f'output/{save_dir}/log_{idx}')
    
    return (popt,perr)  