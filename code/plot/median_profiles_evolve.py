from matplotlib import cm
from matplotlib import pyplot as plt 
from matplotlib.colors import BoundaryNorm
plt.style.use('code/style.mplstyle')
import numpy as np
import argparse
from func import *
import os

# The gravitational softening length (after multiplied with 2.8):
# TNG300: 4kpc
# MTNG: 10.3kpc

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--sim',     default=None,type=str)
parser.add_argument('--DM',      default=None,type=str) # [ Hydro / DM ]
parser.add_argument('--Nboots',  default=1024,type=int)
parser.add_argument('--mass_cut',default=None,type=float)
args = parser.parse_args()

print('')
print('>>> Plot physical density profile evolution <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

str_mass_cut = float_to_str(args.mass_cut)

stats_dir = 'result/bootstrap_stats_phys/'
profile_dir = 'result/bootstrap_phys/'
if args.sim == 'TNG300':
    sim_dir = f'TNG300/sim_205_1250_{args.DM}/'
    snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
    soften_length = 4    # [kpc]
elif args.sim == 'MTNG':
    sim_dir = f'MTNG/{args.DM}-Arepo/MTNG-L500-4320-A/'
    snaps = [264, 237, 214, 179, 151, 129]
    soften_length = 10.3 # [KPC]
    


# Setup the plot
fig, axs = plt.subplots(2, 1, figsize = (5, 10), sharex = True, dpi=500)
    
if args.DM == 'Hydro':
    cmap_name = 'plasma'
elif args.DM == 'DM':
    cmap_name = 'viridis'

num_z = len(snaps)
z_i, z_f = load_z(f'{stats_dir}/{sim_dir}/Nboots_{args.Nboots}/', snaps)
cmap = plt.get_cmap(cmap_name, num_z)
bound = np.linspace(z_f, z_i+0.001, num_z) 
norm = BoundaryNorm(bound, cmap.N)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  ax=axs, orientation='horizontal', spacing='proportional', ticks=bound)
cb.set_label('z')

# Redshift data
for isnap, snap in enumerate(snaps):
        
    # Load data
    data = np.load(f'{stats_dir}/{sim_dir}/Nboots_{args.Nboots}/snap_{snap}_Rsp_stats.npy', allow_pickle = True).item()
    z = data['z']
    scale_factor = 1/(1+z)
    h = 0.6774
    
    mass_bins = data['mass_bins']
    median_idx = data['median_idx_in_boots']
    
    # select median data
    try: 
        mass_cut_idx = np.where(mass_bins == args.mass_cut)[0][0]
        
        profile = np.load(f'{profile_dir}/{sim_dir}/snap_{snap}/Nboots_{args.Nboots}/data/'+
                f'/mass_cut_{str_mass_cut}/boots_{median_idx[mass_cut_idx]}.npy',
                allow_pickle=True).item()
        
        phy_Rsp = data['final_results'][mass_cut_idx, 0, 1] # [kpc]
        phy_R200 = profile['R200_median']                   # [kpc]
        scale_Rsp =  phy_Rsp / phy_R200                     # [dimensionless]
        del data
        
        # Compute softening length in R200
        scale_Rsoft = soften_length / phy_R200              # [dimensionless]

        # density profile
        ax0 = axs[0]
        ax0.set_title(f'{args.sim}-{args.DM} mass cut: 10^{10+args.mass_cut} [$M_\\odot$]')
        # ax0.errorbar(profile['radius'], profile['rho'], yerr = profile['rho_err'].T, color=cmap(norm(np.round(z, 3))),
        #             label=f'z={np.round(z, 3)}') # dimensionless radius
        ax0.plot(profile['fitted_radius'], profile['fitted_rho'], color=cmap(norm(np.round(z, 3))))  # dimensionless radius
        ax0.axvline(x=scale_Rsp, color=cmap(norm(np.round(z, 3))),
                    linestyle='--', linewidth=1)                              # dimensionless radius
        ax0.axvline(x=scale_Rsoft, color=cmap(norm(np.round(z, 3))),
                    linestyle='-', linewidth=1)                               # dimensionless radius
        
        # gradient profile
        ax1 = axs[1]
        # ax1.errorbar(profile['radius'], profile['slope'], yerr = profile['slope_err'].T, # color='r',
        #             label=f'z={np.round(z, 3)}') # dimensionless radius
        ax1.plot(profile['fitted_radius'], profile['fitted_slope'], color=cmap(norm(np.round(z, 3))))  # dimensionless radius
        ax1.axvline(x=scale_Rsp, color=cmap(norm(np.round(z, 3))),
                    linestyle='--', linewidth=1)                              # dimensionless radius
        ax1.axvline(x=scale_Rsoft, color=cmap(norm(np.round(z, 3))),
                    linestyle='-', linewidth=1)                               # dimensionless radius
        
    except IndexError:
        pass 
# general settings
axs[0].set_xticks([]) 
axs[0].set_xscale('log')
axs[1].set_xscale('log')
axs[1].set_xlabel(r"Radius r/$R_{200}$")

axs[0].set_yscale('log')
axs[0].set_ylabel(r"$\rho$/$\rho_c$")
axs[1].set_ylabel(r'd log $\rho$ / d log r')
axs[1].set_ylim(-5, -0.5)

save_dir = f'result/bootstrap_plot_phys/full_{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/{args.sim}-{args.DM}_mass_cut_{str_mass_cut}_profiles')
plt.close