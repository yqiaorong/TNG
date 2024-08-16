import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feat_idx',default=0,type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

print('')
print(f'>>> Plot Rsp feats vs mass <<<')
print('\nInput arguments:')
for key, val in vars(args).items():
	print('{:16} {}'.format(key, val))
print('')

feats = ['Rsp', 'depth', 'width_dimless', 'width_phys']
feat_idx = args.feat_idx

root_dir = 'result/bootstrap_stats/'

# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)

##############################################################################################
### Plot TNG300 ### 
##############################################################################################

TNG300_dir = f'{root_dir}/sim_205_1250/'
TNG300_list = os.listdir(TNG300_dir)

TNG300_cmap = plt.get_cmap('winter', len(TNG300_list))
TNG300_snaps = [99, 78, 67, 50, 40, 33, 21, 17, 13, 8]
TNG300_mass_cuts = [10**11, 10**11.5, 10**12, 10**12.5, 10**13, 10**13.5]

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
    
##############################################################################################
### Plot MTNG-DM ### 
##############################################################################################

MTNG_DM_dir = f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/output'
MTNG_DM_list = os.listdir(MTNG_DM_dir)

MTNG_DM_cmap = plt.get_cmap('autumn', len(TNG300_list))
MTNG_DM_snaps = [264, 237, 214, 179, 151, 129]
MTNG_DM_mass_cuts = [10**13.5, 10**14, 10**14.5, 10**15, 10**15.5]

for isnap, snap in enumerate(MTNG_DM_snaps):
    # Load data
    data = np.load(MTNG_DM_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    z = data['z']
    data = data['final_results']
    num_cut = data.shape[0]
    
    factors = [1000, 1, 1, 1000]
    data[:, feat_idx, :] = data[:, feat_idx, :]*factors[feat_idx] # convert Rsp in [Mpc] to [kpc]
    
    axs.plot(MTNG_DM_mass_cuts[:num_cut], data[:, feat_idx, 1], 
             color=MTNG_DM_cmap(isnap / len(TNG300_snaps)), label=f'z = {np.round(z, 1)}')
    axs.fill_between(MTNG_DM_mass_cuts[:num_cut], data[:, feat_idx, 0], data[:, feat_idx, 2], 
                     color=MTNG_DM_cmap(isnap / len(TNG300_snaps)), alpha=0.2)

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
    axs.set_ylabel(r"Wdith [kpc]")
    axs.set_yscale('log')
axs.legend()

# Save the plot
save_dir = f'result/bootstrap_plot/full/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/dm_{feats[feat_idx]}_vs_mass')
plt.close()