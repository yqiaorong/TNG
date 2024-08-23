import os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse
from func import *

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feat_idx',default=0,   type=int) 
# Feature index [Rsp = 0, depth = 1, width_dimless = 2, width_phys = 3]
parser.add_argument('--variable',default=None,type=str) # [mass / redshift]
args = parser.parse_args()

feat_idx = args.feat_idx
variable = args.variable
if variable is None:
    variable = input("Please provide the variable [M/z]: ")

print('')
print(f'>>> Curve fit as a func of {args.variable} <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Save dir
save_dir = f'result/fitting_partial/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# Run the previous script to get the compiled data
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
os.system(f'python3 code/plot/plot_tgt_mass.py --feat_idx {feat_idx}')

# Load the compiled data
data = np.load(f'data/equ_data_{feats[feat_idx]}.npy', allow_pickle=True).item()



if variable == 'M':
    x = data['mass_x1']
    x = np.log10(x) # put into log space
    constraints = np.round(data['z_x2'], 1)
    xlabel = 'log(Mass)'
elif variable == 'z':
    x = data['z_x2']
    constraints = data['mass_x1']
    xlabel = 'z'
    
uniq_constraints = np.unique(constraints)

y, ymax, ymin = data['y'], data['y_max'], data['y_min']  
if feat_idx == 0  or feat_idx == 3: # [Rsp / width_phys]
    y, ymax, ymin = np.log10(y), np.log10(ymax), np.log10(ymin)

# Print the results setup
if feat_idx == 0:
    output = 'log(R_sp)'
elif feat_idx == 1:
    output = 'depth'
elif feat_idx == 2:
    output = 'width'
elif feat_idx == 3:
    output = 'log(width)'
    
if variable == 'M':
    Input = 'log(M)'
else:
    Input = 'z'
    
# Set up the plot
cmap = plt.get_cmap('autumn', len(uniq_constraints))
fig, axs = plt.subplots(1, len(uniq_constraints), figsize=(15, 2), dpi=500)
ylabels = [r'log($R_{sp}$ [kpc])', 'depth', r'width [$R_{200}$]', 'log(width [kpc])']

# Fit as a func of mass / redshift
for icons, cons in enumerate(uniq_constraints):
    mask = (constraints == cons)
    num_data_points = x[mask].shape[0]
    # print(f'num data points: {num_data_points}')
    
    if num_data_points > 3:

        # fitting
        params, _ = curve_fit(partial_func, x[mask], y[mask], 
                              sigma=(ymax[mask]-ymin[mask])/2, 
                              absolute_sigma=True)
        tot_params[icons] = params
        print(f'equ: {output} ({Input}) = {params[0]:.3f} ({Input})^2 + {params[1]:.3f} ({Input}) + {params[2]:.3f}')
        for key, param in zip(['a:', 'b:', 'c:'], params):
            print('{:2} {}'.format(key, param))
        print('')
        
        # theory
        x_theory = np.linspace(min(x[mask]), max(x[mask]), num=100, endpoint=True)
        y_theory = partial_func(x_theory, *params)
        
        # plot every constraint one by one on one plot
        axs[icons].errorbar(x[mask], y[mask], yerr=[np.abs(ymin[mask]-y[mask]), 
                                                    np.abs(ymax[mask]-y[mask])],
                            fmt='.', c=cmap(icons/len(uniq_constraints)))
        axs[icons].plot(x_theory, y_theory, c='k')
        
        # Set up the labels
        axs[icons].set_xlabel(xlabel)
        if icons == 0:
            axs[icons].set_ylabel(ylabels[feat_idx])
            
plt.savefig(f'{save_dir}/fit_{feats[feat_idx]}_vs_{variable}')
plt.close()

# tot_params = tot_params[~np.all(tot_params == 0, axis=1)] # Remove zeros 
# print('a     b     c')
# print(tot_params)



# Do the one fitting
print('final fitting: ')
params, pcov = curve_fit(partial_func, x, y, sigma=(ymax-ymin)/2, absolute_sigma=True)
print(f'equ: {output} ({Input}) = {params[0]:.3f} ({Input})^2 + {params[1]:.3f} ({Input}) + {params[2]:.3f}')
for key, param in zip(['a:', 'b:', 'c:'], params):
    print('{:2} {}'.format(key, param))
print('')

# theory
x_theory = np.linspace(min(x), max(x), num=100, endpoint=True)
y_theory = partial_func(x_theory, *params)

# plot every constraint together on one plot
fig, ax = plt.subplots(1, 1, dpi=500)
for icons, cons in enumerate(uniq_constraints):
    mask = (constraints == cons)
    
    sorted_idx = np.argsort(x[mask])
    ax.errorbar(x[mask][sorted_idx], y[mask][sorted_idx], 
                yerr=[np.abs(ymin[mask][sorted_idx]-y[mask][sorted_idx]), 
                    np.abs(ymax[mask][sorted_idx]-y[mask][sorted_idx])], 
                fmt='.', c=cmap(icons/len(uniq_constraints)))
    # ax.fill_between(x[mask][sorted_idx], ymax[mask][sorted_idx], ymin[mask][sorted_idx], 
    #                 color=cmap(icons/len(uniq_constraints)), alpha=0.2)
ax.plot(x_theory, y_theory, c='k')

ax.set_xlabel(xlabel)
ax.set_ylabel(ylabels[feat_idx])
plt.savefig(f'{save_dir}/fit_all_{feats[feat_idx]}_vs_{variable}')
plt.close()