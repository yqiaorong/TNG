"""The script uses data from both TNG300 and MTNG to plot how the splashback 
features change with mass. The results are saved in ./result/bootstrap_plot/full/"""

import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',default='',   type=str)
parser.add_argument('--Nboots',  default=None,type=int)
parser.add_argument('--feat_idx',default=0,   type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print('>>> Plot Rsp feats vs mass <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Initial settings
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
feat_idx = args.feat_idx
print(feats[feat_idx])

if args.DM == '':
    cmap_name = 'autumn'
elif args.DM == '_DM':
    cmap_name = 'winter'
    
    

root_dir = 'result/bootstrap_stats/'

# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)

# Saved data
saved_x1, saved_x2, saved_y, saved_y_max, saved_y_min = [], [], [], [], []

##############################################################################################
### Plot TNG300 ### 
##############################################################################################

TNG300_dir = f'{root_dir}/TNG300/sim_205_1250{args.DM}/Nboots_{args.Nboots}/'
TNG300_list = os.listdir(TNG300_dir)
    
TNG300_cmap = plt.get_cmap(cmap_name, len(TNG300_list))
TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
TNG300_mass_cuts = [10**11, 10**11.5, 10**12, 10**12.5, 10**13, 10**13.5, 10**14, 10**14.5]

for isnap, snap in enumerate(TNG300_snaps):
    
    # Load data
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    z = data['z']
    data = data['final_results']
    num_cut = data.shape[0]

    axs.plot(TNG300_mass_cuts[:num_cut], data[:, feat_idx, 1], 
            color=TNG300_cmap(isnap / len(TNG300_snaps)), label=f'z = {np.round(z, 1)}')
    axs.fill_between(TNG300_mass_cuts[:num_cut], data[:, feat_idx, 0], data[:, feat_idx, 2], 
                    color=TNG300_cmap(isnap / len(TNG300_snaps)), alpha=0.2)

    # Append to saved data
    saved_x1.append(TNG300_mass_cuts[:num_cut])
    saved_x2.append([z]*len(TNG300_mass_cuts[:num_cut]))
    saved_y.append(data[:, feat_idx, 1])
    saved_y_min.append(data[:, feat_idx, 0])
    saved_y_max.append(data[:, feat_idx, 2])

#############################################################################################
# Plot MTNG #
#############################################################################################

# if args.DM == '':
#     MTNG_DM_dir = f'{root_dir}/Hydro-Arepo/MTNG-L500-4320-A/output/Nboots_{args.Nboots}/'
# elif args.DM == '_DM':
#     MTNG_DM_dir = f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/output/Nboots_{args.Nboots}/'
# MTNG_DM_list = os.listdir(MTNG_DM_dir)

# MTNG_DM_cmap = plt.get_cmap(cmap_name, len(TNG300_list))
# MTNG_DM_snaps = [264, 237, 214, 179, 151, 129]
# MTNG_DM_mass_cuts = [10**13, 10**13.5, 10**14, 10**14.5, 10**15, 10**15.5]

# for isnap, snap in enumerate(MTNG_DM_snaps):
    
#     # Load data
#     data = np.load(MTNG_DM_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
#     z = data['z']
#     data = data['final_results']
#     num_cut = data.shape[0]
    
#     factors = [1000, 1, 1, 1000]
#     data[:, feat_idx, :] = data[:, feat_idx, :]*factors[feat_idx] # convert Rsp in [Mpc] to [kpc]
    
#     axs.plot(MTNG_DM_mass_cuts[:num_cut], data[:, feat_idx, 1], 
#              color=MTNG_DM_cmap(isnap / len(TNG300_snaps)), label=f'z = {np.round(z, 1)}')
#     axs.fill_between(MTNG_DM_mass_cuts[:num_cut], data[:, feat_idx, 0], data[:, feat_idx, 2], 
#                      color=MTNG_DM_cmap(isnap / len(TNG300_snaps)), alpha=0.2)

#     # Append to saved data
#     saved_x1.append(MTNG_DM_mass_cuts[:num_cut])
#     saved_x2.append([z]*len(MTNG_DM_mass_cuts[:num_cut]))
#     saved_y.append(data[:, feat_idx, 1])
#     saved_y_min.append(data[:, feat_idx, 0])
#     saved_y_max.append(data[:, feat_idx, 2])

# Final edit
axs.set_xscale('log')
axs.set_xlabel('Mass [$M_\\odot$/h]')
if feat_idx == 0:
    axs.set_ylabel(r"$R_{sp}$ [kpc]")
    axs.set_yscale('log')
elif feat_idx == 1:
    axs.set_ylabel("depth")
elif feat_idx == 2:
    axs.set_ylabel(r'width [$R_{200}$]')
elif feat_idx == 3:
    axs.set_ylabel(r"wdith [kpc]")
    axs.set_yscale('log')
axs.legend()

# Save the plot
save_dir = f'result/bootstrap_plot/full{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/dm_{feats[feat_idx]}_vs_mass')
plt.close()



# Flatten the data
saved_x1 = np.array(np.concatenate(saved_x1).tolist())
saved_x2 = np.array(np.concatenate(saved_x2).tolist())
saved_y  = np.array(np.concatenate(saved_y).tolist())
saved_y_min = np.array(np.concatenate(saved_y_min).tolist())
saved_y_max = np.array(np.concatenate(saved_y_max).tolist())

# Save the data
np.save(f'data/equ_data_{feats[feat_idx]}_Nboots{args.Nboots}{args.DM}', 
        {'mass_x1': saved_x1, 'z_x2': saved_x2,
         'y': saved_y, 'y_min': saved_y_min, 'y_max': saved_y_max})
print('data saved.')



# # Curve fit
# from scipy.optimize import curve_fit

# def model_func(x, a, b, c, d, e, f):
#     logM, z = x
#     return a*logM**2 + b*logM*z + c*z**2 + d*logM + e*z + f

# x = np.vstack((np.log10(saved_x1), saved_x2))
# if feat_idx == 0 or 3:
#     param, pcov = curve_fit(model_func, x, np.log10(saved_y), 
#                             sigma=(np.log10(saved_y_max)-np.log10(saved_y_min))/2, 
#                             absolute_sigma=True)
#     print(f'func: log SP (logM, z) = {param[0]:.3f} (log M)^2 + {param[1]:.3f} (log M) z '+
#       f'+ {param[2]:.3f} z^2 + {param[3]:.3f} logM + {param[4]:.3f} z + {param[5]:.3f}')
# else:
#     param, pcov = curve_fit(model_func, x, saved_y, 
#                             sigma=(saved_y_max-saved_y_min)/2, 
#                             absolute_sigma=True)
#     print(f'equ: SP (logM, z) = {param[0]:.3f} (log M)^2 + {param[1]:.3f} (log M) z '+
#         f'+ {param[2]:.3f} z^2 + {param[3]:.3f} logM + {param[4]:.3f} z + {param[5]:.3f}')
# print(param)