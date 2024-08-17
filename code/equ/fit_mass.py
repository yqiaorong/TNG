import os
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse

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
    input = 'log(M)'
else:
    input = 'z'
    
# Set up the plot
cmap = plt.get_cmap('autumn', len(uniq_constraints))
fig, axs = plt.subplots(1, len(uniq_constraints), figsize=(15, 2), dpi=500)
ylabels = [r'log($R_{sp}$ [kpc])', 'depth', r'width [$R_{200}$]', 'log(width [kpc])']

# Fit as a func of mass / redshift
def model_func(x, a, b, c):
    return a*x**2 + b*x + c

for icons, cons in enumerate(uniq_constraints):
    mask = (constraints == cons)
    num_data_points = x[mask].shape[0]
    print(f'num data points: {num_data_points}')
    
    if num_data_points > 3:

        # fitting
        params, pcov = curve_fit(model_func, x[mask], y[mask], 
                                sigma=(ymax[mask]-ymin[mask])/2, 
                                absolute_sigma=True)
        
        print(f'equ: {output} ({input}) = {params[0]:.3f} ({input})^2 + {params[1]:.3f} ({input}) + {params[2]:.3f}')
        for key, param in zip(['a:', 'b:', 'c:'], params):
            print('{:2} {}'.format(key, param))
        print('')
        
        # theory
        x_theory = np.linspace(min(x[mask]), max(x[mask]), num=100, endpoint=True)
        y_theory = model_func(x_theory, *params)
        
        # plot
        axs[icons].errorbar(x[mask], y[mask], yerr=[np.abs(ymin[mask]-y[mask]), 
                                                    np.abs(ymax[mask]-y[mask])],
                            fmt='.', c=cmap(icons/len(uniq_constraints)))
        axs[icons].plot(x_theory, y_theory, c='k')
        # Set up the labels
        axs[icons].set_xlabel(xlabel)
        if icons == 0:
            axs[icons].set_ylabel(ylabels[feat_idx])


# Save the plot
save_dir = f'result/fitting_plot/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/fit_{feats[feat_idx]}_vs_{variable}')
plt.close()