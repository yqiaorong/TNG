"""The script uses data from both TNG300 and MTNG to plot how the splashback 
features width change with accretion rate width."""

import os
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
import seaborn as sns
from scipy.stats import pearsonr
plt.style.use('code/style.mplstyle')
from func import *
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',    default=None,type=str) # [ Hydro /DM ]
parser.add_argument('--Nboots',default=1024,type=int)
parser.add_argument('--width', default=None,type=str) # [ std / percentile ]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats (width) vs accretion rate width ({args.width}) <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

# Initial settings
feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
feat_idx = 2
print(feats[feat_idx])
    
    

root_dir = 'result/bootstrap_stats_phys/'

# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)

##############################################################################################
# Load z
##############################################################################################

# Load TNG300
TNG300_dir = f'{root_dir}/TNG300/sim_205_1250_{args.DM}/Nboots_{args.Nboots}/'
TNG300_snaps = [99, 67, 40, 25, 13, 8]
TNG300_z_i, TNG300_z_f = load_z(TNG300_dir, TNG300_snaps)

# Load MTNG
MTNG_dir = f'{root_dir}/MTNG/{args.DM}-Arepo/MTNG-L500-4320-A/Nboots_{args.Nboots}/'
MTNG_snaps = [264, 214, 151]
MTNG_z_i, MTNG_z_f = load_z(MTNG_dir, MTNG_snaps)

z_i, z_f = max(TNG300_z_i, MTNG_z_i), min(TNG300_z_f, MTNG_z_f)
print(z_i, z_f)
num_z = len(TNG300_snaps)

##############################################################################################
# Set up the plot
##############################################################################################

fig, axs = plt.subplots(1, 1, figsize=(6, 5), dpi=500)

if args.DM == 'Hydro':
    cmap_name = 'plasma'
elif args.DM == 'DM':
    cmap_name = 'viridis'

cmap = plt.get_cmap(cmap_name, num_z)
bound = np.linspace(z_f, z_i+0.01, num_z) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label('z')

##############################################################################################
### Load accretion rate ### 
##############################################################################################

TNG300_mass_cuts, TNG300_z, _, TNG300_accret_width = load_accret(f'result/accretion_rate_plot/TNG300/sim_205_1250_{args.DM}/TNG300_{args.DM}_accret_stats.npy', width=args.width)

MTNG_mass_cuts, MTNG_z, _, MTNG_accret_width = load_accret(f'result/accretion_rate_plot/MTNG/{args.DM}-Arepo/MTNG_{args.DM}_accret_stats.npy', width=args.width)

##############################################################################################
### Plot ### 
##############################################################################################
tot_x, tot_y = [], []

x, y = plot_data(TNG300_dir, TNG300_snaps, [TNG300_z, TNG300_mass_cuts, TNG300_accret_width], axs, cmap, norm)
tot_x = np.concatenate((tot_x, x))
tot_y = np.concatenate((tot_y, y))

x, y = plot_data(MTNG_dir, MTNG_snaps, [MTNG_z, MTNG_mass_cuts, MTNG_accret_width], axs, cmap, norm)
tot_x = np.concatenate((tot_x, x))
tot_y = np.concatenate((tot_y, y))

##############################################################################################
### Correlations ### 
##############################################################################################

# sns.set(style="whitegrid")
# r_value, p_value = pearsonr(tot_x, tot_y)
# reg_color = cmap(norm(np.round(4, 3)))
# axs = sns.regplot(x=tot_x, y=tot_y, ci=95, scatter=False, line_kws={"color": reg_color})
# plt.title(f'r = {r_value:.2f}, p = {p_value:.3f}', loc='center')

# Final edit
axs.set_xlabel(f'Accretion rate width {args.width}')
axs.set_ylabel("Splashback feature width")

# Save the plot
save_dir = f'result/bootstrap_plot_phys/full_{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/{args.DM}_sp_width_vs_accret_width_{args.width}')
plt.close()