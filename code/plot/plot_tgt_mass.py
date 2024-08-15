import os
import numpy as np
from matplotlib import pyplot as plt 
plt.style.use('code/style.mplstyle')
import argparse

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feat_idx',default=0,type=int) # Feature index [Rsp = 0, depth = 1, width = 2]
args = parser.parse_args()

feats = feats = ['Rsp', 'depth', 'width']
feat_idx = args.feat_idx

# Load stats results
root_dir = 'result/bootstrap_stats/'
TNG300_dir = f'{root_dir}/sim_205_1250/'
TNG300_list = os.listdir(TNG300_dir)

MTNG_DM_dir = f'{root_dir}/DM-Arepo/MTNG-L500-4320-A/output'
MTNG_DM_list = os.listdir(MTNG_DM_dir)

# Set up the plot
fig, axs = plt.subplots(1, 1, dpi=500)



### Plot TNG300 ### 

cmap = plt.get_cmap('winter', len(TNG300_list))
TNG300_snaps = [8, 13, 17, 21, 25, 33, 40, 50, 67, 78, 99]
TNG300_mass_cuts = [10**11, 10**11.5, 10**12, 10**12.5, 10**13, 10**13.5]

for isnap, snap in enumerate(TNG300_snaps):
    # Load data
    data = np.load(TNG300_dir+f'/snap_{snap}_Rsp_stats.npy', allow_pickle=True).item()
    z = data['z']
    data = data['final_results']
    num_cut = data.shape[0]
    
    axs.plot(TNG300_mass_cuts[:num_cut], data[:, feat_idx, 1], 
             color=cmap(isnap / len(TNG300_snaps)), label=f'z = {np.round(z, 1)}')
    axs.fill_between(TNG300_mass_cuts[:num_cut], data[:, feat_idx, 0], data[:, feat_idx, 2], 
                     color=cmap(isnap / len(TNG300_snaps)), alpha=0.2)

# Final edit
axs.set_xlabel('Mass [$M_\\odot$/h]')
if args.feat_idx == 0:
    axs.set_ylabel(r"$R_{sp}$ [kpc]")
else:
    axs.set_ylabel(feats[feat_idx])
axs.set_xscale('log')
if args.feat_idx == 0:
    axs.set_yscale('log')
axs.legend()

# Save the plot
save_dir = f'result/bootstrap_plot/full/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
plt.savefig(f'{save_dir}/dm_{feats[feat_idx]}_vs_mass')
plt.close()