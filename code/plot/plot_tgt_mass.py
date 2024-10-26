"""The script uses data from both TNG300 and MTNG to plot how the splashback 
features change with mass. The results are saved in ./result/bootstrap_plot/full/"""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
from func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',      default=None,type=str) # [ Hydro /DM ]
parser.add_argument('--Nboots',  default=1024,type=int)
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
    
    

root_dir = 'result/bootstrap_stats_phys/'

# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)

##############################################################################################
# Load z
##############################################################################################

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{args.DM}/Nboots_{args.Nboots}/'
TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]

TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/{args.DM}-Arepo/MTNG-L500-4320-A/Nboots_{args.Nboots}/'
MTNG_snaps = [264, 237, 214, 179, 151, 129]

MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)
num_z = len(TNG300_snaps)

##############################################################################################
# Set up the plot
##############################################################################################

fig, axs = plt.subplots(1, 1, figsize = (6, 5), dpi=500)

if args.DM == 'Hydro':
    cmap_name = 'plasma'
elif args.DM == 'DM':
    cmap_name = 'viridis'

cmap = plt.get_cmap(cmap_name, num_z)
bound = np.linspace(z_f, z_i+0.001, num_z) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label('z')


##############################################################################################
### Plot ### 
##############################################################################################

# Plot TNG300
axs = plot_stats(TNG300_dir, TNG300_snaps, axs, cmap, norm, feat_idx)
# Plot MTNG #
axs = plot_stats(MTNG_dir, MTNG_snaps, axs, cmap, norm, feat_idx)

# Final edit
axs.set_xscale('log')
axs.set_xlabel('Mass [$M_\\odot$]')
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



# Save the plot
save_dir = f'result/bootstrap_plot_phys/full_{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/{args.DM}_{feats[feat_idx]}_vs_mass')
plt.close()