""""This script plots depth as a function of accretion rate.
    Hydro simulation only."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm


plt.style.use('code/style.mplstyle')
from func import *

print('')
print('>>> Plot mass vs accretion rate <<<')
print('')

root_dir = 'result/bootstrap_stats_phys/'

# ============================================================================================
# Load z
# ============================================================================================

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_Hydro/Nboots_1024/'
TNG300_snaps = [99, 67, 40, 25, 13, 8]
TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

all_z = load_all_z(TNG300_dir, TNG300_snaps)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/Hydro-Arepo/MTNG-L500-4320-A/Nboots_1024/'
MTNG_snaps = [264, 214, 151]
MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)
num_z = len(TNG300_snaps)

# ============================================================================================
# Set up the plot
# ============================================================================================

fig, axs = plt.subplots(1, 1, figsize = (4, 4), dpi=500, sharex=True, constrained_layout=True)
cmap = plt.get_cmap('managua', num_z)
bound = all_z
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
# Reduce colormap ticks sf
from matplotlib.ticker import FuncFormatter
def custom_format(x, pos):
    return f'{x:.1f}'  
cb.ax.xaxis.set_major_formatter(FuncFormatter(custom_format)) 
cb.set_label('z')

# ============================================================================================
# Load accretion rate
# ============================================================================================

def customize(mass, z, accret, accret_err):
    num_cut, num_z = len(mass), len(z)

    accret = accret.reshape(-1)
    accret_err = accret_err.reshape(-1)/2
    
    mass = np.array([10**(mass+10)]*num_z).reshape(-1)    
    z = np.repeat(z, num_cut)
    
    # remove zeros
    non_zero_idx = np.where(accret != 0)[0]
    mass, z, accret, accret_err = mass[non_zero_idx], z[non_zero_idx], accret[non_zero_idx], accret_err[non_zero_idx]
    return mass, z, accret, accret_err

# TNG300
TNG300_mass_cuts, TNG300_z, TNG300_accret_med = load_accret(f'result/accretion_rate_plot/TNG300/sim_205_1250_Hydro/TNG300_Hydro_accret_stats.npy')
_, _, _, TNG300_accret_width = load_accret(f'result/accretion_rate_plot/TNG300/sim_205_1250_Hydro/TNG300_Hydro_accret_stats.npy', 
                                                                 width='percentile')
TNG300_mass_cuts, TNG300_z, TNG300_accret_med, TNG300_accret_err = customize(TNG300_mass_cuts, TNG300_z, TNG300_accret_med, TNG300_accret_width)

# MTNG
MTNG_mass_cuts, MTNG_z, MTNG_accret_med = load_accret(f'result/accretion_rate_plot/MTNG/Hydro-Arepo/MTNG_Hydro_accret_stats.npy')
_, _, _, MTNG_accret_width = load_accret(f'result/accretion_rate_plot/MTNG/Hydro-Arepo/MTNG_Hydro_accret_stats.npy', 
                                                           width='percentile')
MTNG_mass_cuts, MTNG_z, MTNG_accret_med, MTNG_accret_err = customize(MTNG_mass_cuts, MTNG_z, MTNG_accret_med, MTNG_accret_width)

# ============================================================================================
# Plot 
# ============================================================================================

# Plot masses vs accretion rate
for z in np.unique(TNG300_z):
    iz = np.where(TNG300_z == z)[0]
    axs.errorbar(TNG300_mass_cuts[iz], TNG300_accret_med[iz], 
                 yerr = MTNG_accret_err[iz]/2,
                 color=cmap(norm(np.round(z, 3))))
for z in np.unique(MTNG_z):
    iz = np.where(MTNG_z == z)[0]
    axs.errorbar(MTNG_mass_cuts[iz], MTNG_accret_med[iz], 
                 yerr = MTNG_accret_err[iz]/2,
                 color=cmap(norm(np.round(z, 3))), ls='--')
axs.set_xscale('log')
axs.set_xlabel(r'$M_{200m}. [M_\odot]$')
axs.set_ylabel(r'$\Gamma$')

# ============================================================================================
# Curve fit
# ============================================================================================

tot_mass_cuts = np.concatenate((TNG300_mass_cuts, MTNG_mass_cuts))
tot_z         = np.concatenate((TNG300_z, MTNG_z))
tot_accret    = np.concatenate((TNG300_accret_med, MTNG_accret_med))
tot_accret_err = np.concatenate((TNG300_accret_err, MTNG_accret_err))

from scipy.optimize import curve_fit

def fit_func(Inputs, 
             # a, b, 
             c, d
             ):
    """Params:
        Inputs: (x, redshifts)
    """
    x, z = Inputs
    return c*x*z**d # + a*x**b

popt, pcov = curve_fit(fit_func, (np.log10(tot_mass_cuts), tot_z), tot_accret, p0=[1]*2, maxfev=10000)
fit_accret = fit_func((np.log10(tot_mass_cuts), tot_z), *popt)

# Calculate the reduced chi-square
red_chi2  = np.sum(((tot_accret - fit_accret)/ tot_accret_err)**2) / (len(tot_accret) - len(popt))
# Plot
axs.scatter(tot_mass_cuts, fit_accret, c='black', marker='x', s=20,
            label=r'$\chi^2_{\nu}$ = '+f'{red_chi2:.4f}')
axs.set_xscale('log')
axs.legend()

# ============================================================================================
# Confidence interval of popt
# ============================================================================================

perr = np.sqrt(np.diag(pcov))

dof = max(1, len(tot_accret)-len(popt)) 

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

# ============================================================================================
# Save the plot
# ============================================================================================

save_dir = f'result/paper_plots/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

plt.savefig(f'{save_dir}/fig10.png')
plt.close()