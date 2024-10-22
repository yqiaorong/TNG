"""The script uses data from both TNG300 and MTNG to plot how the splashback 
features change with redshift. The results are saved in ./result/bootstrap_plot/full/"""

import os
import numpy as np
from func import *
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--z_or_a',  default='z', type=str)
parser.add_argument('--DM',      default=None,type=str)
parser.add_argument('--Nboots',  default=1024,type=int)
parser.add_argument('--feat_idx',default=0,   type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs {args.z_or_a} <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Initial settings
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
feat_idx = args.feat_idx
print(feats[feat_idx])
    
    

root_dir = 'result/bootstrap_stats_phys/'

##############################################################################################
### Loadd data ###
##############################################################################################

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{args.DM}/Nboots_{args.Nboots}/'
TNG300_list = os.listdir(TNG300_dir)

TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]

TNG300_z, TNG300_bins, TNG300_data = load_data(TNG300_dir, TNG300_snaps, args)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/{args.DM}-Arepo/MTNG-L500-4320-A/Nboots_{args.Nboots}/'
MTNG_list = os.listdir(MTNG_dir)

MTNG_snaps = [264, 237, 214, 179, 151, 129]

MTNG_z, MTNG_bins, MTNG_data = load_data(MTNG_dir, MTNG_snaps, args)

# Concatenate all data
z = np.concatenate((TNG300_z, MTNG_z))
bins = np.concatenate((TNG300_bins, MTNG_bins))
bins = np.log10(bins)
data = np.concatenate((TNG300_data, MTNG_data), axis=0)

##############################################################################################
# Set up the plot
##############################################################################################

fig, axs = plt.subplots(1, 1, figsize = (6, 5), dpi=500)

if args.DM == 'Hydro':
    cmap_name = 'plasma'
elif args.DM == 'DM':
    cmap_name = 'viridis'
    
uniq_bins = np.unique(bins)
num_bins = len(uniq_bins)
min_bin, max_bin = np.min(bins), np.max(bins)

cmap = plt.get_cmap(cmap_name, num_bins)
bound = np.linspace(min_bin, max_bin+0.1, num_bins) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label(r'$\log_{10}(\mathrm{Mass})\ [\log_{10}(M_{\odot})]$')

##############################################################################################
# Plot
##############################################################################################

# Iterate over mass bins
for b in uniq_bins:
    # print('mass cut: ', b)
    idx = np.where(bins == b)
    
    current_z = z[idx]
    current_data = data[idx]
    
    # Plot SORTED data
    sorted_idx = np.argsort(current_z)
    axs.plot(current_z[sorted_idx], current_data[:, feat_idx, 1][sorted_idx], color=cmap(norm(b)), alpha=0.2)
    # axs.fill_between(current_z[sorted_idx], 
    #                  current_data[:, feat_idx, 0][sorted_idx], current_data[:, feat_idx, 2][sorted_idx],
    #                  color=cmap(norm(b)), alpha=0.1)
    axs.errorbar(current_z, current_data[:, feat_idx, 1], 
                yerr=[current_data[:, feat_idx, 1]-current_data[:, feat_idx, 0], 
                      current_data[:, feat_idx, 2]-current_data[:, feat_idx, 1]],
                fmt='.', color=cmap(norm(b))
                    )
    
# Final edit
axs.set_xlabel(args.z_or_a)
if feat_idx == 0:
    axs.set_ylabel(r"$R_{sp}$ [kpc]")
    axs.set_yscale('log')
elif feat_idx == 1:
    axs.set_ylabel("depth")
elif feat_idx == 2:
    axs.set_ylabel(r'width [$R_{200}$]')
elif feat_idx == 3:
    axs.set_ylabel(r"Wdith [kpc]")
    axs.set_yscale('log')



# Save the plot
save_dir = f'result/bootstrap_plot_phys/full_{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'{args.DM}_{feats[feat_idx]}_vs_{args.z_or_a}'))
plt.close()