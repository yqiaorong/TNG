"""The plots are saved in ./result/fitting_full/"""

import os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse
# from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',      default='',  type=str)
parser.add_argument('--Nboots',  default=1024,type=int)
parser.add_argument('--feat_idx',default=0,   type=int) 
parser.add_argument('--func_idx',default=0,   type=int) 
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
if args.DM == '':
    os.system(f'python3 code/plot/plot_tgt_mass.py --feat_idx {feat_idx} --Nboots {args.Nboots}')
elif args.DM == '_DM':
    os.system(f'python3 code/plot/plot_tgt_mass.py --DM {args.DM} --feat_idx {feat_idx} --Nboots {args.Nboots}')

### Load the compiled data ###
data = np.load(f'data/equ_data_{feats[feat_idx]}_Nboots{args.Nboots}{args.DM}.npy', allow_pickle=True).item()

logM, z = np.log10(data['mass_x1']), data['z_x2']
uniq_logM, uniq_z = np.unique(logM), np.unique(z)

y, ymax, ymin = data['y'], data['y_max'], data['y_min']
if feat_idx == 0  or feat_idx == 3: # [Rsp / width_phys]
    y, ymax, ymin = np.log10(y), np.log10(ymax), np.log10(ymin)

# Define functions
def func0(x, 
          a, b, 
          h, i, 
          r, s):
    logM, z = x
    part1 = a*logM + b*logM**2
    part2 = h*z    + i*z**2
    part3 = r*logM*z + s
    return part1 + part2 + part3

# def func1(x, a, b, c, d, e, f):
#     logM, z = x
#     return a*logM**2 + b*logM + c + d*z**2 + e*z + f

def func2(x, 
          a, b, c,
          h, i, j,
          p, q, r, s):
    logM, z = x
    part1 = a*logM + b*logM**2 + c*logM**3
    part2 = h*z    + i*z**2    + j*z**3
    part3 = p*logM**2*z + q*logM*z**2 + r*logM*z + s
    return  part1 + part2 + part3

def func3(x, 
          a, b, c, d, 
          h, i, j, k, 
          m, n, o, 
          p, q, r, s):
    logM, z = x
    part1 = a*logM + b*logM**2 + c*logM**3 + d*logM**4
    part2 = h*z    + i*z**2    + j*z**3    + k*z**4
    part3 = m*logM**3*z + n*logM**2*z**2 + o*logM*z**3
    part4 = p*logM**2*z + q*logM*z**2    + r*logM*z + s
    return  part1 + part2 + part3 + part4

def func4(x, 
          a, b, c, d, e,
          h, i, j, k, l,
          m, n, o, 
          p, q, r, s,
          c1, c2, c3, c4):
    logM, z = x
    part1 = a*logM + b*logM**2 + c*logM**3 + d*logM**4 + e*logM**5
    part2 = h*z    + i*z**2    + j*z**3    + k*z**4    + l*z**5
    part3 = m*logM**3 *z + n*logM**2 *z**2 + o*logM    *z**3
    part4 = p*logM**2 *z + q*logM    *z**2 + r*logM    *z    + s
    part5 = c1*logM**4*z + c2*logM**3*z**2 + c3*logM**2*z**3 + c4*logM*z**4
    return  part1 + part2 + part3 + part4 + part5

def func5(x, 
          a, b, c, d, e, f,
          h, i, j, k, l, m1, 
          m, n, o, 
          p, q, r, s,
          c1, c2, c3, c4,
          c5, c6, c7, c8, c9):
    logM, z = x
    part1 = a*logM + b*logM**2 + c*logM**3 + d*logM**4 + e*logM**5 + f*logM**6
    part2 = h*z    + i*z**2    + j*z**3    + k*z**4    + l*z**5    + m1*z**6
    part4 = p *logM**2*z + q *logM   *z**2 + r*logM    *z    + s
    part3 = m *logM**3*z + n *logM**2*z**2 + o*logM    *z**3
    part5 = c1*logM**4*z + c2*logM**3*z**2 + c3*logM**2*z**3 + c4*logM*z**4
    part6 = c5*logM**5*z + c6*logM**4*z**2 + c7*logM**3*z**3 + c8*logM**2*z**4 + c9*logM*z**5 
    return  part1 + part2 + part3 + part4 + part5 + part6

func_list = [func0, # func1, 
             func2, func3, func4, func5]
sel_func = func_list[args.func_idx]

# Full fit
variables = np.vstack((logM, z))
sketchy_factor = 1
mean_err = (ymax-ymin)*sketchy_factor/2
popt, _ = curve_fit(sel_func, variables, y, sigma=mean_err, absolute_sigma=True)
# print(f'func: log SP (logM, z) = {popt[0]:.3f} (log M)^2 + {popt[1]:.3f} (log M) z '+
#     f'+ {popt[2]:.3f} z^2 + {popt[3]:.3f} logM + {popt[4]:.3f} z + {popt[5]:.3f}')
print(popt)

# Validation
def chi_sq(x, y, yerr, func, popt):
    y_pred = func(x, *popt)
    deg_f = (len(y_pred)-len(popt))
    
    chi2_stats = np.sum(((y-y_pred)/yerr)**2) 
    chi2_redu = chi2_stats / deg_f
    return chi2_stats, chi2_redu, deg_f

chi2_stats, chi2_redu, deg_f = chi_sq(variables, y, mean_err, sel_func, popt)

# P value
from scipy.stats import chi2
p_val = 1 - chi2.cdf(chi2_stats, deg_f)
print(f'p value: {p_val}')

########################################################################################
# Plot
########################################################################################

if args.DM == '_DM':
    cmap_name = 'winter'
elif args.DM == '':
    cmap_name = 'autumn'
cmap1 = plt.get_cmap(cmap_name, len(uniq_logM))
cmap2 = plt.get_cmap(cmap_name, len(uniq_z))
fig, axs = plt.subplots(1, 2, figsize=(9, 4), dpi=500)

# As a func of mass
for ifix_z, fix_z in enumerate(uniq_z):
    mask = (z == fix_z)
    logM_th = np.linspace(min(logM[mask]), max(logM[mask]), num=100, endpoint=True)
    # Input to func
    var = (logM_th, np.broadcast_to(fix_z, logM_th.shape))
    # Output of func
    y_th = sel_func(var, *popt)
    
    # Plot the raw data
    axs[0].errorbar(logM[mask], y[mask], yerr=[sketchy_factor*np.abs(ymin[mask]-y[mask]), 
                                               sketchy_factor*np.abs(ymax[mask]-y[mask])],
                    fmt='.', # c=f'C{ifix_z}',
                    c=cmap1(ifix_z/len(uniq_z))
                    )
    # Plot the fitted line
    axs[0].plot(logM_th, y_th, # c=f'C{ifix_z}', 
                c=cmap1(ifix_z/len(uniq_z)), 
                label=f'z = {fix_z:.1f}')
    
for ifix_logM, fix_logM in enumerate(uniq_logM):
    mask = (logM == fix_logM)
    z_th = np.linspace(min(z[mask]), max(z[mask]), num=100, endpoint=True)
    # Input to func
    var = (np.broadcast_to(fix_logM, z_th.shape), z_th)
    # Output of func
    y_th = sel_func(var, *popt)
    # Plot the raw data
    axs[1].errorbar(z[mask], y[mask], yerr=[sketchy_factor*np.abs(ymin[mask]-y[mask]), 
                                            sketchy_factor*np.abs(ymax[mask]-y[mask])],
                    fmt='.', # c=f'C{ifix_logM}',
                    c=cmap2(ifix_logM/len(uniq_logM))
                    )
    # Plot the fitted line
    axs[1].plot(z_th, y_th, # c=f'C{ifix_logM}', 
                c=cmap2(ifix_logM/len(uniq_logM)), 
                label=rf'mass = $10^{{{fix_logM:.1f}}}$'+'$M_\\odot$/h')
    plt.text(0.73, 0.1, f'$\\chi^2_{{{{red}}}} = {chi2_redu:.2f}$', transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top')
    plt.text(0.73, 0.15, f'p val = {p_val:.3f}', transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top')

# Label
axs[0].set_xlabel('log(Mass [$M_\\odot$/h])')
axs[1].set_xlabel('z')
if args.func_idx == 2 and args.feat_idx == 2:
    axs[0].set_ylim(top=3)
    axs[1].set_ylim(top=3)
if args.feat_idx == 2:
    axs[0].set_ylim(bottom=0)
    axs[1].set_ylim(bottom=0)
ylabels = [r'log($R_{sp}$ [kpc])', 'depth', r'width [$R_{200}$]', 'log(width [kpc])']
axs[0].set_ylabel(ylabels[feat_idx])
axs[0].legend()
axs[1].legend()

# Save 
if args.DM == '':
    save_type = 'Hydro'
elif args.DM == '_DM':
    save_type = 'DM'
save_dir = f'result/fitting_full/{save_type}-sim_Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/{save_type}_fit_{feats[feat_idx]}_func{args.func_idx}')
plt.close()