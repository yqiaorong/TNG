import os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse
# from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feat_idx',default=0, type=int) 
parser.add_argument('--func_idx',default=0, type=int) 
# Feature index [Rsp = 0, depth = 1, width_dimless = 2, width_phys = 3]
args = parser.parse_args()

print('')
print(f'>>> Curve fit as a func of mass and z <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))

# Run the previous script to get the compiled data
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
feat_idx = args.feat_idx
os.system(f'python3 code/plot/plot_tgt_mass.py --feat_idx {feat_idx}')

### Load the compiled data ###
data = np.load(f'data/equ_data_{feats[feat_idx]}.npy', allow_pickle=True).item()

logM, z = np.log10(data['mass_x1']), data['z_x2']
uniq_logM, uniq_z = np.unique( logM), np.unique(z)

y, ymax, ymin = data['y'], data['y_max'], data['y_min'] 
if feat_idx == 0  or feat_idx == 3: # [Rsp / width_phys]
    y, ymax, ymin = np.log10(y), np.log10(ymax), np.log10(ymin)

# Define functions
def func1(x, a, b, c, d, e, f):
    logM, z = x
    return a*logM**2 + b*logM*z + c*z**2 + d*logM + e*z + f

def func2(x, a, b, c, d, e, f):
    logM, z = x
    return a*logM**2 + b*logM + c + d*z**2 + e*z + f

func_list = [func1, func2]
sel_func = func_list[args.func_idx]

# Full fit
variables = np.vstack((logM, z))
popt, _ = curve_fit(sel_func, variables, y, sigma=(ymax-ymin)/2, absolute_sigma=True)
print(f'func: log SP (logM, z) = {popt[0]:.3f} (log M)^2 + {popt[1]:.3f} (log M) z '+
    f'+ {popt[2]:.3f} z^2 + {popt[3]:.3f} logM + {popt[4]:.3f} z + {popt[5]:.3f}')
print(popt)

# Validation
def chi_sq(x, y, yerr, func, popt):
    y_pred = func(x, *popt)
    return np.sum(((y-y_pred)/yerr)**2) / (len(y_pred)-len(popt))

chi = chi_sq(variables, y, (ymax-ymin)/2, sel_func, popt)

########################################################################################
# Plot
########################################################################################

cmap1 = plt.get_cmap('autumn', len(uniq_logM))
cmap2 = plt.get_cmap('cool', len(uniq_z))
fig, axs = plt.subplots(1, 2, figsize=(6, 3), dpi=500)

# As a func of mass
for ifix_z, fix_z in enumerate(uniq_z):
    mask = (z == fix_z)
    logM_th = np.linspace(min(logM[mask]), max(logM[mask]), num=100, endpoint=True)
    # Input to func
    var = (logM_th, np.broadcast_to(fix_z, logM_th.shape))
    # Output of func
    y_th = sel_func(var, *popt)
    # Plot the raw data
    axs[0].errorbar(logM[mask], y[mask], yerr=[np.abs(ymin[mask]-y[mask]), 
                                               np.abs(ymax[mask]-y[mask])],
                    fmt='.', c=cmap1(ifix_z/len(uniq_z)))
    # Plot the fitted line
    axs[0].plot(logM_th, y_th, c=cmap1(ifix_z/len(uniq_z)), label=f'z = {fix_z:.1f}')
    
for ifix_logM, fix_logM in enumerate(uniq_logM):
    mask = (logM == fix_logM)
    z_th = np.linspace(min(z[mask]), max(z[mask]), num=100, endpoint=True)
    # Input to func
    var = (np.broadcast_to(fix_logM, z_th.shape), z_th)
    # Output of func
    y_th = sel_func(var, *popt)
    # Plot the raw data
    axs[1].errorbar(z[mask], y[mask], yerr=[np.abs(ymin[mask]-y[mask]), 
                                            np.abs(ymax[mask]-y[mask])],
                    fmt='.', c=cmap2(ifix_logM/len(uniq_logM)))
    # Plot the fitted line
    axs[1].plot(z_th, y_th, c=cmap2(ifix_logM/len(uniq_logM)), 
                label=rf'mass = $10^{{{fix_logM:.1f}}}$'+'$M_\\odot$/h')
    plt.text(0.05, 0.95, f'$\\chi^2 = {chi:.2f}$', transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top')

# Label
axs[0].set_xlabel('log(Mass [$M_\\odot$/h])')
axs[1].set_xlabel('z')
ylabels = [r'log($R_{sp}$ [kpc])', 'depth', r'width [$R_{200}$]', 'log(width [kpc])']
axs[0].set_ylabel(ylabels[feat_idx])
axs[0].legend()
axs[1].legend()

# Save 
save_dir = f'result/fitting_full/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/fit_{feats[feat_idx]}_func{args.func_idx}')
plt.close()

