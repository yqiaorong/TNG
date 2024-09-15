"""The script uses data from both TNG300 and MTNG to plot how the splashback 
features change with redshift. The results are saved in ./result/bootstrap_plot/full/"""

import os
import h5py
import numpy as np
import illustris_python as il
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--DM',  default='',type=str)
parser.add_argument('--Nboots',  default=None,type=int)
parser.add_argument('--feat_idx',default=0,   type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs mass <<<')
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

##############################################################################################
### Plot TNG300 ### 
##############################################################################################

TNG300_dir = f'{root_dir}/TNG300/sim_205_1250{args.DM}/Nboots_{args.Nboots}/'
TNG300_list = os.listdir(TNG300_dir)

TNG300_cmap = plt.get_cmap(cmap_name, len(TNG300_list))
TNG300_snaps = [99, 78, 67, 50, 40, 33, 25, 21, 17, 13, 8]
TNG300_mass_cuts = [11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5]

TNG300_z, TNG300_all_data = [], []
for snap in TNG300_snaps: # from low z to high z (present)
        
    # Load data
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    z = data['z']
    TNG300_z.append(z)
    data = data['final_results']
    TNG300_all_data.append(data)

for cut_idx in range(len(TNG300_mass_cuts)-1): # from low cut to high cut
                     
    # Sort all data
    plot_x, plot_y, plot_y_min, plot_y_max = [], [], [], []
    
    for snap_idx in range(len(TNG300_all_data)):
        if TNG300_all_data[snap_idx].shape[0] > cut_idx:
            plot_y.append(TNG300_all_data[snap_idx][cut_idx, feat_idx, 1])
            plot_y_min.append(TNG300_all_data[snap_idx][cut_idx, feat_idx, 0])
            plot_y_max.append(TNG300_all_data[snap_idx][cut_idx, feat_idx, 2])
            
            plot_x.append(TNG300_z[snap_idx])
        else:
            pass
 
    axs.plot(plot_x, plot_y, color=TNG300_cmap(cut_idx / len(TNG300_mass_cuts)),
             label=r'$10^{%.1f}$'%TNG300_mass_cuts[cut_idx]+'~'
             +r'$10^{%.1f}$ '%TNG300_mass_cuts[cut_idx+1]+'$M_\\odot$/h')
    axs.fill_between(plot_x, plot_y_min, plot_y_max, 
                     color=TNG300_cmap(cut_idx / len(TNG300_mass_cuts)), alpha=0.2)
    
##############################################################################################
### Plot MTNG-DM ### 
##############################################################################################

# if args.DM == '':
#     MTNG_DM_dir = f'{root_dir}/Hydro-Arepo/MTNG-L500-4320-A/output/Nboots_{args.Nboots}/'
# elif args.DM == '_DM':
#     MTNG_DM_dir = f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/output/Nboots_{args.Nboots}/'
# MTNG_DM_list = os.listdir(MTNG_DM_dir)

# MTNG_DM_cmap = plt.get_cmap(cmap_name, len(TNG300_list))
# MTNG_DM_snaps = [264, 237, 214, 179, 151, 129]
# MTNG_DM_mass_cuts = [13, 13.5, 14, 14.5, 15]

# MTNG_DM_z, MTNG_DM_all_data = [], []
# for snap in MTNG_DM_snaps: # from high z to low z (present)
        
#     # Load data
#     data = np.load(MTNG_DM_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
#     z = data['z']
#     MTNG_DM_z.append(z)
#     data = data['final_results']
#     MTNG_DM_all_data.append(data)

# factors = [1000, 1, 1, 1000]
# for cut_idx in range(len(MTNG_DM_mass_cuts)-1): # from low cut to high cut
    
#     # Because of difference in unit, Rsp should be multipled by 1000 from [Mpc] to [kpc]
#     factor = factors[feat_idx]
                     
#     # Sort all data
#     plot_x, plot_y, plot_y_min, plot_y_max = [], [], [], []
    
#     for snap_idx in range(len(MTNG_DM_all_data)):
#         if MTNG_DM_all_data[snap_idx].shape[0] > cut_idx:
#             plot_y.append(factor*MTNG_DM_all_data[snap_idx][cut_idx, feat_idx, 1])
#             plot_y_min.append(factor*MTNG_DM_all_data[snap_idx][cut_idx, feat_idx, 0])
#             plot_y_max.append(factor*MTNG_DM_all_data[snap_idx][cut_idx, feat_idx, 2])
            
#             plot_x.append(MTNG_DM_z[snap_idx])
#         else:
#             pass

#     axs.plot(plot_x, plot_y, color=MTNG_DM_cmap(cut_idx / len(MTNG_DM_mass_cuts)),
#              label=r'$10^{%.1f}$'%MTNG_DM_mass_cuts[cut_idx]+'~'
#              +r'$10^{%.1f}$ '%MTNG_DM_mass_cuts[cut_idx+1]+'$M_\\odot$/h')
#     axs.fill_between(plot_x, plot_y_min, plot_y_max, 
#                      color=MTNG_DM_cmap(cut_idx / len(MTNG_DM_mass_cuts)), alpha=0.2)

# Final edit
axs.set_xlabel('z')
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
axs.legend()

# Save the plot
save_dir = f'result/bootstrap_plot/full{args.DM}/Nboots_{args.Nboots}/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(os.path.join(save_dir, f'dm_{feats[feat_idx]}_vs_redshift'))
plt.close()